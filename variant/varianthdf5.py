from collections import defaultdict
from gzip import open
from pickle import dump, load
from pprint import pprint

from tables import (Filters, Float32Col, HDF5ExtError, Int32Col, IsDescription,
                    StringCol, open_file)

from .variant import (get_ann, get_genotype,
                      get_variant_start_and_end_positions, get_variant_type)

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
        self.id_to_chrom_dict_file_path = '{}.id_to_chrom_dict.pickle.gz'.format(
            self.variant_file_path)

        # Data
        self.variant_hdf5 = None
        self.id_to_chrom_dict = {}

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
        Initialize self.variant_hdf5 & self.id_to_chrom_dict.
        :param reset: bool; re-make data instead of reading from files
        :return: None
        """

        if not reset:

            try:
                print('Reading variant HDF5 ...')
                self.variant_hdf5 = open_file(
                    self.variant_hdf5_file_path, mode='r')

                print('Reading ID-to-chromosome dict ...')
                self._read_id_to_chrom_dict()

            except (FileNotFoundError, HDF5ExtError) as e:
                print('\tFailed: {}.'.format(e))
                reset = True

        if reset:
            print('Resetting ...')

            if self.variant_hdf5:
                self.variant_hdf5.close()
                print('Closed varinat HDF5.')

            print('Making variant HDF5 ...')
            self._make_variant_hdf5()

            print('Reading variant HDF5 ...')
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

                    # Parse .VCF row
                    chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample = line.split(
                        '\t')[:10]

                    start, end = get_variant_start_and_end_positions(
                        int(pos), ref, alt)

                    variant_type = get_variant_type(ref, alt)

                    # pathogenicity = None

                    effect, impact, gene_name = get_ann(
                        ['effect', 'impact', 'gene_name'], info=info)

                    gt = get_genotype(format_, sample)

                    if chrom not in chrom_table_to_row_dict:  # Make table
                        chrom_table = variant_hdf5.create_table(
                            '/',
                            'chromosome_{}_variants'.format(chrom),
                            description=self._VariantDescription,
                            expectedrows=chrom_n_rows[chrom], )
                        print('\t\tMaking chromosome {} variant table ...'.
                              format(chrom))
                        chrom_table_to_row_dict[chrom] = chrom_table.row

                    # Write variant
                    cursor = chrom_table_to_row_dict[chrom]

                    cursor['CHROM'] = chrom
                    cursor['ID'] = id_
                    cursor['REF'] = ref
                    cursor['ALT'] = alt
                    cursor['QUAL'] = qual
                    cursor['start'] = start
                    cursor['end'] = end
                    cursor['variant_type'] = variant_type
                    # cursor['pathogenicity'] = pathogenicity
                    cursor['effect'] = effect
                    cursor['impact'] = impact
                    cursor['gene_name'] = gene_name
                    cursor['GT'] = gt

                    cursor.append()

                    if id_ != '.':
                        self.id_to_chrom_dict[id_] = chrom

                print('Flushing tables and making column indices ...')
                for chrom in chrom_table_to_row_dict:
                    print('\tchromosome {} variant table ...'.format(chrom))
                    chrom_table = variant_hdf5.get_node(
                        '/', 'chromosome_{}_variants'.format(chrom))
                    chrom_table.flush()

                    for col in [
                            'CHROM',
                            'ID',
                            'REF',
                            'ALT',
                            'QUAL',
                            'start',
                            'end',
                            'variant_type',
                            # 'pathogenicity'
                            'effect',
                            'impact',
                            'gene_name',
                            'GT',
                    ]:
                        chrom_table.cols._f_col(col).create_csindex()

                self.variant_hdf5 = variant_hdf5
                print(self.variant_hdf5)

                print('Writing ID-to-chromosome dict ...')
                self._write_id_to_chrom_dict()

    class _VariantDescription(IsDescription):
        """
        Describe VariantHDF5 table columns.
        """

        # TODO: Match with VCF specification
        CHROM = StringCol(8)
        ID = StringCol(16)
        REF = StringCol(256)
        ALT = StringCol(256)
        QUAL = Float32Col()
        # POS
        start = Int32Col()
        end = Int32Col()
        # REF & ALT
        variant_type = StringCol(8)
        # INFO
        pathogenicity = StringCol(8)
        # INFO ANN
        effect = StringCol(64)
        impact = StringCol(8)
        gene_name = StringCol(16)
        # FORMAT & SAMPLE
        GT = StringCol(16)

    def _write_id_to_chrom_dict(self):
        """
        Write ID-to-chromosome dict to file.
        :return: None
        """

        with open(self.id_to_chrom_dict_file_path, 'wb') as f:
            dump(self.id_to_chrom_dict, f)

    def _read_id_to_chrom_dict(self):
        """
        Read ID-to-chromosome dict from file.
        :return: None
        """

        with open(self.id_to_chrom_dict_file_path, 'rb') as f:
            self.id_to_chrom_dict = load(f)

    def get_variants_by_id(self, id_):
        """
        Get variants by id.
        :param id_: str; ID
        :return: dict; {}
        """

        chrom = self.id_to_chrom_dict.get(id_)

        if chrom:  # Variant found
            chrom_table = self.variant_hdf5.get_node(
                '/', 'chromosome_{}_variants'.format(chrom))
            return self._read_where(chrom_table, "ID == b'{}'".format(id_))

        else:  # Variant not found
            return None

    def get_variants_by_region(self, chrom, start, end):
        """
        Get variants by region.
        :param hdf5_table: HDF5 Table
        :param chrom: str; chromosome
        :param start: int; start position
        :param end: int; end position
        :return: list; of dict; (n_results); [{column: value, ...}, ...]
        """

        chrom_table = self.variant_hdf5.get_node(
            '/', 'chromosome_{}_variants'.format(chrom))
        return self._read_where(
            chrom_table, '({} <= start) & (end <= {})'.format(start, end))

    def _read_where(self, hdf5_table, query):
        """
        Do hdf5_table.read_where(query) and map the results to column names.
        :param hdf5_table: HDF5 Table
        :param query: str; query
        :return: list; of dict; (n_results); [{column: value, ...}, ...]
        """

        columns = hdf5_table.colnames

        print('Reading {} where: {} ...'.format(hdf5_table.name, query))
        dicts = []
        for row in hdf5_table.read_where(query):

            dict_ = {}
            for c, v in zip(columns, row):
                try:
                    dict_[c] = v.decode()
                except AttributeError:
                    dict_[c] = v
            dicts.append(dict_)

        return dicts
