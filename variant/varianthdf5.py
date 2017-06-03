from collections import defaultdict
from gzip import open
from pickle import dump, load

from tables import (Filters, Float32Col, HDF5ExtError, Int32Col, IsDescription,
                    StringCol, open_file)

from .hdf5.hdf5.hdf5 import read_where_and_map_column_names
from .vcf import (get_vcf_info, get_vcf_info_ann, get_vcf_sample_format,
                  update_vcf_variant_dict)


class VariantHDF5:
    """
    Data structure storing variants in .HDF5.
    """

    def __init__(self, vcf_file_path, reset=False):
        """
        Construct VariantHDF5.
        :param vcf_file_path: str; .VCF file path
        :param reset: bool; re-make data instead of reading from files
        :return: None
        """

        # File paths
        self.vcf_file_path = vcf_file_path
        self.variant_hdf5_file_path = '{}.hdf5'.format(self.vcf_file_path)
        self.id_to_chrom_dict_file_path = '{}.id_to_chrom_dict.pickle.gz'.format(
            self.vcf_file_path)
        self.gene_to_chrom_dict_file_path = '{}.gene_to_chrom_dict.pickle.gz'.format(
            self.vcf_file_path)

        # Data
        self.variant_hdf5 = None
        self.id_to_chrom_dict = {}
        self.gene_to_chrom_dict = {}

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
        Initialize self.variant_hdf5 & self.id_to_chrom_dict &
        self.gene_to_chrom_dict.
        :param reset: bool; re-make data instead of reading from files
        :return: None
        """

        if not reset:

            try:
                print('Reading ...')

                print('\tVariant HDF5 ...')
                self.variant_hdf5 = open_file(
                    self.variant_hdf5_file_path, mode='r')

                print('\tID-to-chromosome dict ...')
                self._read_id_to_chrom_dict()

                print('\tgene-to-chromosome dict ...')
                self._read_gene_to_chrom_dict()

            except (OSError, FileNotFoundError, HDF5ExtError) as e:
                print('\tFailed: {}.'.format(e))
                reset = True

        if reset:
            print('Resetting ...')

            if self.variant_hdf5:
                self.variant_hdf5.close()
                print('Closed variant HDF5.')

            self._make_variant_hdf5()

            print('Reading variant HDF5 ...')
            self.variant_hdf5 = open_file(
                self.variant_hdf5_file_path, mode='r')

    def _make_variant_hdf5(self):
        """
        Make variant .HDF5.
        :return: None
        """

        with open(self.vcf_file_path, 'rt') as f:

            print('Getting data-start position ...')
            data_start_position = None
            line = f.readline()
            while line.startswith('#'):
                data_start_position = f.tell()
                line = f.readline()

            print('Counting variants per chromosome ...')
            chrom_n_rows = defaultdict(lambda: 0)
            chrom = None
            while line:

                a_chrom = line.split('\t')[0]
                if a_chrom != chrom:
                    print('\t@ {} ...'.format(a_chrom))
                    chrom = a_chrom
                chrom_n_rows[a_chrom] += 1

                line = f.readline()

            print('Making variant HDF5 ...')
            with open_file(
                    self.variant_hdf5_file_path,
                    mode='w',
                    filters=Filters(
                        complevel=1, complib='blosc')) as variant_hdf5:

                # Cursors for chromosome tables
                chrom_table_to_row_dict = {}

                f.seek(data_start_position)
                for i, line in enumerate(f):
                    if i % 1000000 == 0:
                        print('\t@ {} ...'.format(i + 1))

                    # Parse .VCF row
                    chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample = line.split(
                        '\t')

                    if chrom not in chrom_table_to_row_dict:  # Make table

                        chrom_table = variant_hdf5.create_table(
                            '/',
                            'chromosome_{}_variants'.format(chrom),
                            description=self._VariantDescription,
                            expectedrows=chrom_n_rows[chrom])
                        print('\t\tMaking {} table ...'.format(
                            chrom_table.name))

                        chrom_table_to_row_dict[chrom] = chrom_table.row

                    # Write variant
                    cursor = chrom_table_to_row_dict[chrom]

                    cursor['CHROM'] = chrom
                    cursor['POS'] = pos
                    cursor['ID'] = id_
                    cursor['REF'] = ref
                    cursor['ALT'] = alt
                    cursor['QUAL'] = qual
                    caf = get_vcf_info('CAF', info=info)
                    if caf:
                        cursor['CAF'] = caf
                    clnsig = get_vcf_info('CLNSIG', info=info)
                    if clnsig:
                        cursor['CLNSIG'] = clnsig
                    clndbn = get_vcf_info('CLNDBN', info=info)
                    if clndbn:
                        try:
                            cursor['CLNDBN'] = clndbn
                        except TypeError:
                            print('\tCLNDBN error with {}'.format(clndbn))
                    cursor['effect'] = get_vcf_info_ann('effect', info=info)[0]
                    cursor['impact'] = get_vcf_info_ann('impact', info=info)[0]
                    gene_name = get_vcf_info_ann('gene_name', info=info)[0]
                    cursor['gene_name'] = gene_name
                    cursor['GT'] = get_vcf_sample_format(
                        'GT', format_=format_, sample=sample)

                    cursor.append()

                    # Update *-to-chromosome dict
                    if id_ != '.':
                        self.id_to_chrom_dict[id_] = chrom

                    self.gene_to_chrom_dict[gene_name] = chrom

                print('\tFlushing tables and making column indices ...')
                for chrom in chrom_table_to_row_dict:
                    chrom_table = variant_hdf5.get_node(
                        '/', 'chromosome_{}_variants'.format(chrom))
                    print('\t\t{} table ...'.format(chrom_table.name))
                    chrom_table.flush()

                    for col in [
                            'CHROM',
                            'POS',
                            'ID',
                            'REF',
                            'ALT',
                            'QUAL',
                            'CAF',
                            'CLNSIG',
                            'CLNDBN',
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

                print('Writing gene-to-chromosome dict ...')
                self._write_gene_to_chrom_dict()

    class _VariantDescription(IsDescription):
        """
        Describe VariantHDF5 table columns.
        """

        # TODO: Optimize
        CHROM = StringCol(8)
        POS = Int32Col()
        ID = StringCol(16)
        REF = StringCol(256)
        ALT = StringCol(256)
        QUAL = Float32Col()
        # INFO
        CAF = StringCol(16)
        CLNSIG = StringCol(8)
        CLNDBN = StringCol(256)
        # INFO ANN
        effect = StringCol(32)
        impact = StringCol(32)
        gene_name = StringCol(32)
        # FORMAT & sample
        GT = StringCol(8)

    def _read_id_to_chrom_dict(self):
        """
        Read ID-to-chromosome dict from file.
        :return: None
        """

        with open(self.id_to_chrom_dict_file_path, 'rb') as f:
            self.id_to_chrom_dict = load(f)

    def _write_id_to_chrom_dict(self):
        """
        Write ID-to-chromosome dict to file.
        :return: None
        """

        with open(self.id_to_chrom_dict_file_path, 'wb') as f:
            dump(self.id_to_chrom_dict, f)

    def _read_gene_to_chrom_dict(self):
        """
        Read gene-to-chromosome dict from file.
        :return: None
        """

        with open(self.gene_to_chrom_dict_file_path, 'rb') as f:
            self.gene_to_chrom_dict = load(f)

    def _write_gene_to_chrom_dict(self):
        """
        Write gene-to-chromosome dict to file.
        :return: None
        """

        with open(self.gene_to_chrom_dict_file_path, 'wb') as f:
            dump(self.gene_to_chrom_dict, f)

    def get_variants_by_id(self, id_):
        """
        Get variants by id.
        :param id_: str; ID
        :return: dict; variant dict
        """

        chrom = self.id_to_chrom_dict.get(id_)

        if chrom:  # Variant found

            chrom_table = self.variant_hdf5.get_node(
                '/', 'chromosome_{}_variants'.format(chrom))

            variant_dicts = read_where_and_map_column_names(
                chrom_table, "ID == b'{}'".format(id_))

            if 1 < len(variant_dicts):
                raise ValueError('Found 1 < variants with id {}.'.format(id_))

            variant_dict = variant_dicts.pop()
            self._make_variant_dict_consistent(variant_dict)
            update_vcf_variant_dict(variant_dict)

            return variant_dict

    def get_variants_by_gene(self, gene):
        """
        Get variants by gene.
        :param gene: str; HGNC gene name
        :return: list; of variant dict
        """

        chrom = self.gene_to_chrom_dict.get(gene)

        if chrom:  # Variant found

            chrom_table = self.variant_hdf5.get_node(
                '/', 'chromosome_{}_variants'.format(chrom))

            variant_dicts = read_where_and_map_column_names(
                chrom_table, "gene_name == b'{}'".format(gene))

            for d in variant_dicts:
                self._make_variant_dict_consistent(d)
                update_vcf_variant_dict(d)

            return variant_dicts

    def get_variants_by_region(self, chrom, start, end):
        """
        Get variants by region (chrom:start-end).
        :param hdf5_table: HDF5 Table
        :param chrom: str; chromosome
        :param start: int; start position
        :param end: int; end position
        :return: list; of variant dict
        """

        chrom_table = self.variant_hdf5.get_node(
            '/', 'chromosome_{}_variants'.format(chrom))

        variant_dicts = read_where_and_map_column_names(
            chrom_table, '({} <= POS) & (POS <= {})'.format(start, end))

        for d in variant_dicts:
            self._make_variant_dict_consistent(d)
            update_vcf_variant_dict(d)

        return variant_dicts

    def _make_variant_dict_consistent(
            self,
            variant_dict,
            ann_fields=('effect', 'impact', 'gene_name'),
            sample_fields=('GT', )):
        """
        Update .VCF variant dict in place.
        :param dict; variant dict
        :param ann_fields: tuple; variant_dict['ANN'] = {...}
        :param sample_fields: tuple; variant_dict['sample'][0] = {...}
        :return: None
        """

        # Remove field with value: None | ''
        for k in list(variant_dict.keys()):
            if not variant_dict[k]:
                del variant_dict[k]

        variant_dict['ANN'] = {
            0: {field: variant_dict.pop(field)
                for field in ann_fields}
        }

        variant_dict['sample'] = {
            0: {field: variant_dict.pop(field)
                for field in sample_fields}
        }
