from collections import defaultdict
from gzip import open
from pickle import dump, load
from pprint import pprint

from tables import (Filters, Float32Col, HDF5ExtError, Int32Col, IsDescription,
                    StringCol, open_file)

from .variant import get_variant_start_and_end_positions

HDF5_COMPRESSOR = 'blosc'
HDF5_COMPRESSION_LEVEL = 1


class VariantHDF5:
    """
    Data structure storing variants in .HDF5.
    """

    def __init__(self, variant_file_path, reset=False):
        """
        Construct VariantHDF5.
        :param reset: bool; re-make data instead of reading from files
        :return: None
        """

        # File paths
        self.variant_file_path = variant_file_path
        self.variant_hdf5_file_path = '{}.hdf5'.format(self.variant_file_path)
        self.rsis_to_chrom_dict_file_path = '{}.rsid_to_chrom_dict.pickle.gz'.format(
            self.variant_file_path)

        # Data
        self.variant_hdf5 = None
        self.rsid_to_chrom_dict = {}

        self._load_data(reset=reset)

    def __del__(self):
        """
        Destruct VariantHDF5.
        :return: None
        """

        if self.variant_hdf5:
            self.variant_hdf5.close()
            print('Closed variant HDF5.')

    def _load_data(self, reset=False):
        """
        Initialize self.variant_hdf5 & self.rsid_to_chrom_dict.
        :param reset: bool; re-make data instead of reading from files
        :return: None
        """

        if not reset:

            try:
                print('Reading variant HDF5 ...')
                self.variant_hdf5 = open_file(
                    self.variant_hdf5_file_path, mode='r')

                print('Reading RSID-to-chromosome dict ...')
                self._read_rsid_to_chrom_dict()

            except (FileNotFoundError, HDF5ExtError) as e:
                print('\tFailed ({}).'.format(e))
                reset = True

        if reset:
            print('Resetting ...')

            if self.variant_hdf5:
                self.variant_hdf5.close()

            self._make_variant_hdf5()

            self.variant_hdf5 = open_file(
                self.variant_hdf5_file_path, mode='r')

    def _make_variant_hdf5(self):
        """
        Make .HDF5 storing variants.
        :return: None
        """

        with open(self.variant_file_path, 'rt') as f:

            print('Getting data start position ...')
            data_start_position = None
            line = f.readline()
            while line.startswith('#'):
                data_start_position = f.tell()
                line = f.readline()

            print('Counting variants in chromosomes ...')
            chrom_n_rows = defaultdict(lambda: 0)
            chrom = None
            while line:
                a_chrom = line.split('\t')[0]
                if a_chrom != chrom:
                    print('\t@ {} ...'.format(a_chrom))
                    chrom = a_chrom
                    chrom_n_rows[a_chrom] += 1
                line = f.readline()
            pprint(chrom_n_rows)

            print('Making variant HDF5 ...')
            with open_file(
                    self.variant_hdf5_file_path,
                    mode='w',
                    filters=Filters(
                        complevel=HDF5_COMPRESSION_LEVEL,
                        complib=HDF5_COMPRESSOR)) as variant_hdf5:

                chrom_table_to_row_dict = {}

                f.seek(data_start_position)
                for i, line in enumerate(f):
                    if i % 1000000 == 0:
                        print('\t@ {} ...'.format(i + 1))

                    chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample = line.split(
                        '\t')[:10]
                    start, end = get_variant_start_and_end_positions(
                        int(pos), ref, alt)

                    variant_type = get_variant_type(ref, alt)
                    gt = get_genotype(format_, sample)
                    effect, impact, gene_name = get_ann(
                        info, ['effect', 'impact', 'gene_name'])
                    pathogenicity =

                    if chrom not in chrom_table_to_row_dict:
                        chrom_table = variant_hdf5.create_table(
                            '/',
                            'chromosome_{}_variants'.format(chrom),
                            description=self._VariantDescription,
                            #                             expectedrows=chrom_n_rows[chrom],
                        )
                        print('\t\tMaking chromosome {} variant table ...'.
                              format(chrom))
                        chrom_table_to_row_dict[chrom] = chrom_table.row

                    cursor = chrom_table_to_row_dict[chrom]
                    cursor['chrom'] = chrom
                    cursor['start'] = start
                    cursor['end'] = end
                    cursor['id_'] = id_
                    cursor['ref'] = ref
                    cursor['alt'] = alt
                    cursor['qual'] = qual
                    cursor['variant_type'] = variant_type
                    cursor['gt'] = gt
                    cursor['effect'] = effect
                    cursor['impact'] = impact
                    cursor['gene_name'] = gene_name
                    cursor['pathogenicity'] = pathogenicity

                    if id_ != '.':
                        self.rsid_to_chrom_dict[id_] = chrom

                    cursor.append()

                print('Flushing tables and making column indices ...')
                for chrom in chrom_table_to_row_dict:
                    print('\tchromosome {} variant table ...'.format(chrom))
                    chrom_table = variant_hdf5.get_node(
                        '/', 'chromosome_{}_variants'.format(chrom))
                    chrom_table.flush()

                    # TODO: Sort like .VCF row
                    for col in ['chrom', 'start', 'end', 'id_', 'ref', 'alt']:
                        print('\t\t Making {} index ...'.format(col))
                        chrom_table.cols._f_col(col).create_csindex()

                self.variant_hdf5 = variant_hdf5
                print(self.variant_hdf5)

                print('Writing RSID-to-chromosome dict ...')
                self._write_rsid_to_chrom_dict()

    class _VariantDescription(IsDescription):
        """
        Describe VariantHDF5 table columns.
        """

        # TODO: Match with VCF specification
        chrom = StringCol(16)
        start = Int32Col()
        end = Int32Col()
        id_ = StringCol(16)
        ref = StringCol(256)
        alt = StringCol(256)
        qual = Float32Col()
        gt = StringCol(16)
        vt = StringCol(8)
        effect = StringCol(64)
        impact = StringCol(8)
        gene_name = StringCol(16)

    def _write_rsid_to_chrom_dict(self):
        """
        Write RSID-to-chromosome dict to file.
        :return: None
        """

        with open(self.rsis_to_chrom_dict_file_path, 'wb') as f:
            dump(self.rsid_to_chrom_dict, f)

    def _read_rsid_to_chrom_dict(self):
        """
        Read RSID-to-chromosome dict from file.
        :return: None
        """

        with open(self.rsis_to_chrom_dict, 'rb') as f:
            self.rsid_to_chrom_dict = load(f)

    def get_variant_by_id(self, id_):
        """
        Search for id_ in variants.
        :param id_: str; RSID
        :return: dict; {}
        """

        chrom = self.rsid_to_chrom_dict.get(id_)

        if chrom:  # Variant found
            chrom_table = self.variant_hdf5.get_node(
                '/', 'chromosome_{}_variants'.format(chrom))
            return [
                row.fetch_all_fields()
                for row in chrom_table.where("id_ == b'{}'".format(id_))
            ]

        else:  # Variant not found
            return None
