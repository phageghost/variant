# ======================================================================================================================
# DataFrame functions
# ======================================================================================================================
import collections
import datetime
import gzip
import pickle
import re

import tables
from numpy import argmin, empty
from pandas import DataFrame, read_csv

import tabix
from __init__ import (CHROMOSOMES, CHROMOSOMES_CHR, MUTATION_EFFECT_RANKING,
                      PATH_CHAIN_GRCH37_TO_GRCH38, PATH_CHAIN_HG19_TO_HG38,
                      PATH_CHROMOSOME_MAP, PATH_CLINVAR, PATH_DBSNP,
                      PATH_GRCH38, PATH_HG38, PICARD, SNPEFF, SNPSIFT)
from compression import bgzip_tabix
from support.file import mark_filename, read_vcf
from support.system import run_command

VCF_FIELD_CASTER = {
    'POS': int,
    'QUAL': float,
    'GT': lambda x: re.split('[|/]', x),
    'AD': lambda x: x.split(','),
    'VQSLOD': float,
    'CLNSIG': lambda x: max([int(s) for s in re.split('[,|]', x)]),
}

ANN_FIELDS = [
    'ALT', 'effect', 'impact', 'gene_name', 'gene_id', 'feature_type',
    'feature_id', 'transcript_biotype', 'rank', 'hgvsc', 'hgvsp',
    'cdna_position', 'cds_position', 'protein_position', 'distance_to_feature',
    'error'
]


# ======================================================================================================================
# DataFrame functions
# ======================================================================================================================
def convert_vcf_to_maf(vcf,
                       ensg_to_entrez,
                       sample_name=None,
                       n_anns=1,
                       maf_filepath=None):
    """
    .VCF ==> .MAF.
    :param vcf: DataFrame or str;
    :param ensg_to_entrez: str; File path to a mapping compression
    (ENSG ID<\t>Entrez ID)
    :param sample_name:
    :param n_anns: int;
    :param maf_filepath: str; filepath to MAF.
    :return: DataFrame; MAF
    """

    if isinstance(vcf, DataFrame):  # DataFrame
        if not sample_name:  # Sample_name must be provided
            raise ValueError('If variant is a DataFrame, provide sample_name.')

    elif isinstance(vcf, str):  # Filepath to a VCF
        # Read VCF and get sample_name
        vcf_dict = read_vcf(vcf)
        # TODO: handle multisamples
        sample_name = vcf_dict['samples'][0]
        vcf = vcf_dict['data']
        assert isinstance(vcf, DataFrame)

    else:
        raise ValueError(
            'variant must be either a DataFrame or filepath to a VCF.')

    # Read mapping as dictionary
    ensg_to_entrez = read_csv(ensg_to_entrez, index_col=0).to_dict()

    maf_header = [
        'Hugo_Symbol',
        'Entrez_Gene_Id',
        'Center',
        'NCBI_Build',
        'Chromosome',
        'Start_Position',
        'End_Position',
        'Strand',
        'Variant_Classification',
        'Variant_Type',
        'Reference_Allele',
        'Tumor_Seq_Allele1',
        'Tumor_Seq_Allele2',
        'dbSNP_RS',
        'dbSNP_Val_Status',
        'Tumor_Sample_Barcode',
        'Matched_Norm_Sample_Barcode',
        'Matched_Norm_Seq_Allele1',
        'Matched_Norm_Seq_Allele2',
        'Tumor_Validation_Allele1',
        'Tumor_Validation_Allele2',
        'Match_Norm_Validation_Allele1',
        'Match_Norm_Validation_Allele2',
        'Verification_Status',
        'Validation_Status',
        'Mutation_Status',
        'Sequencing_Phase',
        'Sequence_Source',
        'Validation_Method',
        'Score',
        'BAM_File',
        'Sequencer',
        'Tumor_Sample_UUID',
        'Matched_Norm_Sample_UUID',
    ]
    maf = DataFrame(index=vcf.index, columns=maf_header)

    print('Iterating through VCF rows ...')
    tmp = empty((vcf.shape[0], 10), dtype=object)
    for i, vcf_row in vcf.iterrows():  # For each VCF row
        if i % 10000 == 0:
            print('\t VCF row {} ...'.format(i))

        chrom, pos, rsid, ref, alt, info = vcf_row.iloc[[0, 1, 2, 3, 4, 7]]

        start, end = _get_variant_start_and_end_positions(pos, ref, alt)

        if (len(ref) - len(alt)) % 3:
            inframe = False
        else:
            inframe = True

        vt = _get_variant_type(ref, alt)

        for i_s in info.split(';'):  # For each INFO
            if i_s.startswith('ANN='):  # ANN
                # Strip 'ANN=' prefix
                i_s = i_s[4:]
                anns = i_s.split(',')

                for a in anns[:n_anns]:  # For each ANN
                    a_s = a.split('|')

                    effect = a_s[1]

                    vc = _get_maf_variant_classification(effect, vt, inframe)

                    gene_name = a_s[3]

                    gene_id = ensg_to_entrez.get(a_s[4], 0)

        tmp[i] = gene_name, gene_id, chrom, start, end, vc, vt, rsid, ref, alt

    maf.ix[:, [
        'Hugo_Symbol',
        'Entrez_Gene_Id',
        'Chromosome',
        'Start_Position',
        'End_Position',
        'Variant_Classification',
        'Variant_Type',
        'dbSNP_RS',
        'Reference_Allele',
        'Tumor_Seq_Allele1',
    ]] = tmp
    maf.ix[:, 'Strand'] = '+'
    maf.ix[:, [
        'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
        'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID'
    ]] = sample_name

    # Save
    if maf_filepath:
        if not maf_filepath.endswith('.maf'):
            maf_filepath += '.maf'
        maf.to_csv(maf_filepath, sep='\t', index=None)

    return maf


def get_variant_type(vcf_data):
    """

    :param vcf_data: DataFrame;
    :return: Series;
    """

    def f(vcf_row):
        ref, alt = vcf_row.iloc[[3, 4]]
        return _get_variant_type(ref, alt)

    s = vcf_data.apply(f, axis=1)
    s.name = 'variant_type'
    return s


def get_variant_start_and_end_positions(vcf_data):
    """

    :param vcf_data: DataFrame;
    :return: list; list of lists which contain
    """

    def f(vcf_row):
        pos, ref, alt = vcf_row.iloc[[1, 3, 4]]
        return _get_variant_start_and_end_positions(pos, ref, alt)

    s = vcf_data.apply(f, axis=1)
    s.name = 'start_end_positions'
    return s


def get_ann(vcf_data, field, n_anns=1):
    """

    :param vcf_data: DataFrame;
    :param field: str;
    :param n_anns: int;
    :return: list; list of lists which contain field
    """

    i = ANN_FIELDS.index(field)

    def f(vcf_row):
        for i_s in vcf_row.iloc[7].split(';'):  # For each INFO

            if i_s.startswith('ANN='):  # ANN
                # Strip 'ANN=' prefix
                i_s = i_s[4:]

                to_return = []
                anns = i_s.split(',')

                for a in anns[:n_anns]:  # For each ANN
                    to_return.append(a.split('|')[i])

                if len(to_return) == 1:
                    return to_return[0]
                else:
                    return to_return

    s = vcf_data.apply(f, axis=1)
    s.name = field
    return s


# TODO: handle multisample
def get_allelic_frequencies(vcf_data, sample_iloc=9, join=True):
    """

    :param vcf_data: DataFrame;
    :param sample_iloc: int;
    :param join: bool;
    :return: list; of lists, or str if join=True, containing allelic
    frequencies for a sample
    """

    def f(vcf_row):
        s_split = vcf_row.iloc[sample_iloc].split(':')
        try:
            dp = int(s_split[2])
            return tuple([
                '{:0.2f}'.format(ad / dp)
                for ad in [int(i) for i in s_split[1].split(',')]
            ])
        except ValueError:
            return ()
        except ZeroDivisionError:
            return ()

    s = vcf_data.apply(f, axis=1)

    if join:
        s = s.apply(lambda t: ','.join(t))
    s.name = 'allelic_frequency'

    return s


# ======================================================================================================================
# Tabix VCF
# ======================================================================================================================
def get_variants_by_tabix(sample_vcf,
                          contig=None,
                          start=None,
                          end=None,
                          query_str=None,
                          reference_vcf=None):
    """

    :param sample_vcf: str or pytabix handler;
    :param contig: str;
    :param start: int;
    :param end: int;
    :param query_str: int;
    :param reference_vcf: str or pytabix handler;
    :return: list; list of dict
    """

    if isinstance(sample_vcf, str):  # Open sample VCF
        sample_vcf = tabix.open(sample_vcf)

    if query_str:
        records = sample_vcf.querys(query_str)
    else:
        records = sample_vcf.query(contig, start, end)

    if reference_vcf and len(
            list(records)
    ) == 0:  # If sample does not have the record, query reference if given

        if isinstance(reference_vcf, str):  # Open reference VCF
            reference_vcf = tabix.open(reference_vcf)

        records = reference_vcf.query(contig, start - 1, end)

    return [parse_variant(r) for r in records]


# ======================================================================================================================
# Parse VCF fields
# ======================================================================================================================
def parse_variant(vcf_row, n_anns=1):
    """

    :param vcf_row: iterable;
    :param n_anns: int;
    :return: dict; won't have the field if the its value is missing
    """

    variant = {
        'CHROM': _cast_vcf_field('CHROM', vcf_row[0]),
        'POS': _cast_vcf_field('POS', vcf_row[1]),
        'REF': _cast_vcf_field('REF', vcf_row[3]),
        'FILTER': _cast_vcf_field('FILTER', vcf_row[6]),
    }
    ref = variant['REF']

    # ID
    rsid = vcf_row[2]
    if rsid and rsid != '.':
        variant['ID'] = _cast_vcf_field('ID', rsid)

    # ALT
    alt = vcf_row[4]
    if alt and alt != '.':
        variant['ALT'] = _cast_vcf_field('ALT', alt)

    # QUAL
    qual = vcf_row[5]
    if qual and qual != '.':
        variant['QUAL'] = _cast_vcf_field('QUAL', qual)

    # Variant type
    if alt:
        vt = _get_variant_type(ref, alt)
        variant['variant_type'] = vt

    # Samples
    variant['samples'] = []
    format_ = vcf_row[8].split(':')
    for i, s in enumerate(vcf_row[9:]):
        s_d = {'sample_id': i + 1}

        # Sample
        for k, v in zip(format_, s.split(':')):
            s_d[k] = _cast_vcf_field(k, v)

        # Genotype
        if 'ALT' in variant:
            ref_alts = [variant['REF']] + variant['ALT'].split(',')
            s_d['genotype'] = [ref_alts[int(gt)] for gt in s_d['GT']]
        else:
            s_d['genotype'] = [variant['REF']] * 2

        # Allelic frequency
        if 'DP' in s_d and int(s_d['DP']):
            s_d['allelic_frequency'] = [
                round(int(ad) / int(s_d['DP']), 3) for ad in s_d['AD']
            ]

        variant['samples'].append(s_d)

    info_split = vcf_row[7].split(';')
    for i_s in info_split:
        if i_s.startswith('ANN='):
            anns = {}
            for i, a in enumerate(i_s.split(',')[:n_anns]):
                a_split = a.split('|')

                anns[i] = {
                    'effect': a_split[1],
                    'putative_impact': a_split[2],
                    'gene_name': a_split[3],
                    'gene_id': a_split[4],
                    'feature_type': a_split[5],
                    'feature_id': a_split[6],
                    'transcript_biotype': a_split[7],
                    'rank': a_split[8],
                    'hgvsc': a_split[9],
                    'hgvsp': a_split[10],
                    'cdna_position': a_split[11],
                    'cds_position': a_split[12],
                    'protein_position': a_split[13],
                    'distance_to_feature': a_split[14],
                }
            variant['ANN'] = anns
        else:
            try:
                k, v = i_s.split('=')
                if v and v != '.':
                    # TODO: decode properly
                    variant[k] = _cast_vcf_field(k, v)
            except ValueError:
                pass
                # print('INFO error: {} (not key=value)'.format(i_s))

    return variant


def _cast_vcf_field(k, v, caster=VCF_FIELD_CASTER):
    """

    :param k: str;
    :param v: str;
    :param caster: dict;
    :return:
    """

    if k in caster:
        return caster[k](v)
    else:
        return v


def _get_variant_type(ref, alt):
    """

    :param ref: str;
    :param alt: str;
    :return: str;
    """

    if len(ref) == len(alt):
        if len(ref) == 1:
            vt = 'SNP'
        elif len(ref) == 2:
            vt = 'DNP'
        elif len(ref) == 3:
            vt = 'TNP'
        else:  # 4 <= len(ref)
            vt = 'ONP'

    elif len(ref) < len(alt):
        vt = 'INS'

    else:  # len(alt) < len(ref)
        vt = 'DEL'

    return vt


def _get_variant_start_and_end_positions(pos, ref, alt):
    """

    :param ref: str;
    :param alt: str;
    :return: (str, str);
    """

    if len(ref) == len(alt):
        s, e = pos, pos + len(alt) - 1

    elif len(ref) < len(alt):
        s, e = pos, pos + 1

    else:  # len(alt) < len(ref)
        s, e = pos + 1, pos + len(ref) - len(alt)

    return s, e


def _get_maf_variant_classification(es, vt, inframe):
    """

    :param es: str; effect or effects concatenated by '&'
    :param vt: str; Variant type
    :param inframe: bool;
    :return: str; MAF variant classification
    """

    es = es.split('&')
    vc = convert_ann_effect_to_maf_variant_classification(
        es[argmin([MUTATION_EFFECT_RANKING.index(e) for e in es])], vt,
        inframe)
    return vc


def convert_ann_effect_to_maf_variant_classification(e, vt, inframe):
    """

    :param e: str; Effect
    :param vt: str; Variant type
    :param inframe: bool;
    :return: str; MAF variant classification
    """

    if e in (
            'transcript_ablation',
            'exon_loss_variant',
            'splice_acceptor_variant',
            'splice_donor_variant', ):
        vc = 'Splice_Site'

    elif e in ('stop_gained', ):
        vc = 'Nonsense_Mutation'

    elif vt == 'INS' and (e == 'frameshift_variant' or (not inframe and e in (
            'protein_protein_contact',
            'protein_altering_variant', ))):
        vc = 'Frame_Shift_Ins'

    elif vt == 'DEL' and (e == 'frameshift_variant' or (not inframe and e in (
            'protein_protein_contact',
            'protein_altering_variant', ))):
        vc = 'Frame_Shift_Del'

    elif e in ('stop_lost', ):
        vc = 'Nonstop_Mutation'

    elif e in (
            'start_lost',
            'initiator_codon_variant', ):
        vc = 'Translation_Start_Site'

    elif vt == 'INS' and inframe and e in (
            'protein_protein_contact',
            'disruptive_inframe_insertion',
            'inframe_insertion',
            'protein_altering_variant', ):
        vc = 'In_Frame_Ins'

    elif vt == 'DEL' and inframe and e in (
            'protein_protein_contact',
            'disruptive_inframe_deletion',
            'inframe_deletion',
            'protein_altering_variant', ):
        vc = 'In_Frame_Del'

    elif e in (
            'transcript_variant',
            'conservative_missense_variant',
            'rare_amino_acid_variant',
            'missense_variant',
            'coding_sequence_variant', ) or (vt not in ('INS', 'DEL') and
                                             e == 'protein_protein_contact'):
        vc = 'Missense_Mutation'

    elif e in (
            'transcript_amplification',
            'splice_region_variant',
            'intragenic_variant',
            'conserved_intron_variant',
            'intron_variant',
            'INTRAGENIC', ):
        vc = 'Intron'

    elif e in (
            'incomplete_terminal_codon_variant',
            'start_retained_variant',
            'stop_retained_variant',
            'synonymous_variant',
            'NMD_transcript_variant', ):
        vc = 'Silent'

    elif e in (
            'exon_variant',
            'mature_miRNA_variant',
            'non_coding_exon_variant',
            'non_coding_transcript_exon_variant',
            'non_coding_transcript_variant',
            'nc_transcript_variant', ):
        vc = 'RNA'

    elif e in (
            '5_prime_UTR_variant',
            '5_prime_UTR_premature_start_codon_gain_variant', ):
        vc = '5\'UTR'

    elif e in ('3_prime_UTR_variant', ):
        vc = '3\'UTR'

    elif e in (
            'TF_binding_site_ablation',
            'TFBS_ablation',
            'TF_binding_site_amplification',
            'TFBS_amplification',
            'TF_binding_site_variant',
            'TFBS_variant',
            'regulatory_region_ablation',
            'regulatory_region_amplification',
            'regulatory_region_variant',
            'regulatory_region',
            'feature_elongation',
            'feature_truncation',
            'conserved_intergenic_variant',
            'intergenic_variant',
            'intergenic_region', ):
        vc = 'IGR'

    elif e in ('upstream_gene_variant', ):
        vc = '5\'Flank'

    elif e in ('downstream_gene_variant', ):
        vc = '3\'Flank'

    elif e in ('sequence_feature', ):
        vc = 'Targeted_Region'

    else:
        print('Unknown: e={}, vt={}, & inframe={}.'.format(e, vt, inframe))
        vc = 'Targeted_Region'

    return vc


def convert_maf_variant_classification_to_mutsig_effect(vc):
    d = {
        '3\'-UTR': 'noncoding',
        '3\'Flank': 'noncoding',
        '3\'Promoter': 'noncoding',
        '3\'UTR': 'noncoding',
        '5\'-Flank': 'noncoding',
        '5\'-UTR': 'noncoding',
        '5\'Flank': 'noncoding',
        '5\'Promoter': 'noncoding',
        '5\'UTR': 'noncoding',
        'De_novo_Start': 'null',
        'De_novo_Start_InFrame': 'null',
        'De_novo_Start_OutOfFrame': 'null',
        'Frame_Shift_Del': 'null',
        'Frame_Shift_Ins': 'null',
        'IGR': 'noncoding',
        'In_Frame_Del': 'null',
        'In_Frame_Ins': 'null',
        'Intron': 'noncoding',
        'Missense': 'nonsilent',
        'Missense_Mutation': 'nonsilent',
        'NCSD': 'noncoding',
        'Non-coding_Transcript': 'noncoding',
        'Nonsense': 'null',
        'Nonsense_Mutation': 'null',
        'Nonstop_Mutation': 'null',
        'Promoter': 'noncoding',
        'RNA': 'noncoding',
        'Read-through': 'null',
        'Silent': 'silent',
        'Splice': 'null',
        'Splice_Region': 'null',
        'Splice_Site': 'null',
        'Splice_Site_DNP': 'null',
        'Splice_Site_Del': 'null',
        'Splice_Site_Ins': 'null',
        'Splice_Site_ONP': 'null',
        'Splice_Site_SNP': 'null',
        'Start_Codon_DNP': 'null',
        'Start_Codon_Del': 'null',
        'Start_Codon_Ins': 'null',
        'Start_Codon_ONP': 'null',
        'Stop_Codon_DNP': 'null',
        'Stop_Codon_Del': 'null',
        'Stop_Codon_Ins': 'null',
        'Synonymous': 'silent',
        'Targeted_Region': 'silent',
        'Translation_Start_Site': 'null',
        'Variant_Classification': 'effect',
        'downstream': 'noncoding',
        'miRNA': 'noncoding',
        'upstream': 'noncoding',
        'upstream;downstream': 'noncoding',
    }

    return d[vc]


VCF_FILE_PATH = './clinvar_20170501.vcf.gz'
VCF_RSID_INDEX_FILE_PATH = './variant_index'

HDF5_FILE_PATH = './variant_table.hdf5'
HDF5_COMPRESSOR = 'blosc'
HDF5_COMPRESSION_LEVEL = 1


class Variants:
    """
    Serve up variants.
    """

    class VariantTable(tables.IsDescription):
        """
        Table description.
        """

        rsid = tables.StringCol(16)

        # TODO: Consider removing
        is_snp = tables.BoolCol()

        contig = tables.StringCol(64)
        start = tables.Int32Col()
        end = tables.Int32Col()

        ref = tables.StringCol(256)
        alt = tables.StringCol(256)

        qual = tables.Float32Col()

        # TODO: Consider removing
        affected_area = tables.StringCol(64)

        # TODO: Add missing fields
        putative_impact = tables.StringCol(8)
        affected_gene = tables.StringCol(16)
        pathogenicity = tables.Int32Col()

        gt = tables.StringCol(16)

    def __init__(self,
                 vcf_file_path=VCF_FILE_PATH,
                 hdf5_file_path=HDF5_FILE_PATH,
                 vcf_rsid_index_file_path=VCF_RSID_INDEX_FILE_PATH,
                 force_rebuild=False):
        """
        Create an object used to interface with variant tables on disk.
        :param vcf_file_path: str;
        :param hdf5_file_path: str;
        :param vcf_rsid_index_file_path: str;
        :param force_rebuild: bool;
        """

        # Attributes
        self.vcf_file_path = vcf_file_path
        # TODO: Refactor
        self.vcf_field_model = {
            'primary_key': {
                'rsid': (2, str),
            },
            'core_fields': {
                'contig': (0, str),
                'start': (1, int),
                'ref': (3, str),
                'alt': (4, str),
                'qual': (5, float),
            },
            'info_position': 7,
            'secondary_field_keys': 8,
            'secondary_field_values': 9,
        }
        self.vcf_rsid_index_file_path = vcf_rsid_index_file_path
        # TODO: Rename
        self.hdf5 = None
        self.rsid_index = {}
        self.variant_group = None

        # Start
        start_time = datetime.datetime.now()

        if not force_rebuild:

            try:
                print('Opening HDF5 {} ...'.format(hdf5_file_path))
                self.hdf5 = tables.open_file(filename=hdf5_file_path, mode='r')
                print('\tDone in {}.'.format(datetime.datetime.now() -
                                             start_time))
                self._load_rsid_index(self.vcf_rsid_index_file_path)

            except (OSError, IOError, tables.HDF5ExtError):
                print('\tFailed.')
                force_rebuild = True
                if self.hdf5:
                    self.hdf5.close()

        if force_rebuild:  # Load variants from VCF

            self._populate_from_vcf(vcf_file_path, hdf5_file_path)
            self.hdf5 = tables.open_file(filename=hdf5_file_path, mode='r')
            self._save_rsid_index(self.vcf_rsid_index_file_path)

        # TODO: Rename 'variants'
        self.variant_group = self.hdf5.get_node('/variants')

    def __del__(self):
        """
        Destructor. Close HDF5.
        :return: None
        """

        if self.hdf5:
            print('Closing HDF5 ...')
            self.hdf5.close()
        print('\tDone.')

    def _save_rsid_index(self, vcf_rsid_index_file_path):
        """
        Saves the index mapping rsids to contigs to a pickle contained in vcf_rsid_index_file_path
        :param vcf_rsid_index_file_path: str;
        :return: None
        """

        start_time = datetime.datetime.now()

        print('Saving rsid index to {}...'.format(vcf_rsid_index_file_path))
        with gzip.open(vcf_rsid_index_file_path, 'wb') as index_file:
            pickle.dump(
                self.rsid_index, index_file, protocol=-1, fix_imports=True)

        print('\tDone in {}'.format(datetime.datetime.now() - start_time))

    def _load_rsid_index(self, vcf_rsid_index_file_path):
        """
        Load the index mapping rsids to contigs into memory from vcf_rsid_index_file_path
        :param vcf_rsid_index_file_path: str;
        :return: None
        """

        start_time = datetime.datetime.now()
        print('Loading rsid index from {}...'.format(vcf_rsid_index_file_path))
        with gzip.open(vcf_rsid_index_file_path, 'rb') as index_file:
            self.rsid_index = pickle.load(index_file, fix_imports=True)
        print('\tDone in {}'.format(datetime.datetime.now() - start_time))

    def _populate_from_vcf(self, vcf_file_path, hdf5_file_path):
        """
        Generate internal data structures from VCF file and cache them to disk for future reference.
        """

        start_time = datetime.datetime.now()
        print('Parsing VCF and building variant HDF5 ...')

        # Open the VCF and iterate over it
        data_start_pos = None
        with gzip.open(vcf_file_path, 'rt') as vcf_file:
            # Skip past comments and header
            line = vcf_file.readline()
            while line.startswith('#'):
                # Remember this position so we can seek back to it
                data_start_pos = vcf_file.tell()
                line = vcf_file.readline()

            # Count rows
            print('Pass 1 of 2: Counting number of rows in each contig ...')
            precount_start_time = datetime.datetime.now()

            rows_per_contig = collections.defaultdict(lambda: 0)
            line = vcf_file.readline()
            while line != '':
                rows_per_contig[convert_chromosome(
                    line.split('\t')[self.vcf_field_model['core_fields'][
                        'contig'][0]])] += 1
                line = vcf_file.readline()

            print('\tDone in {}.'.format(datetime.datetime.now() -
                                         precount_start_time))
            vcf_file.seek(data_start_pos)

            # Initialize table file
            with tables.open_file(
                    filename=hdf5_file_path,
                    mode='w',
                    title='Customer variants',
                    filters=tables.Filters(
                        complevel=HDF5_COMPRESSION_LEVEL,
                        complib=HDF5_COMPRESSOR)) as h5file:
                # Initialize group
                self.variant_group = h5file.create_group(
                    '/', 'variants', 'Client sequence variants')
                # Create a dictionary of cursors (pointer to the currently active row)
                variant_cursors = {}

                # Parse the tuples in each row and stick them into our tables
                print(
                    'Pass 2 of 2: Populating variant table from VCF contents ...'
                )
                for row_num, vcf_row in enumerate(vcf_file):
                    if row_num % 1000000 == 0:
                        print('\tProcessing VCF row {} ...'.format(row_num +
                                                                   1))
                    if vcf_row is not '':
                        split_row = vcf_row.rstrip().split('\t')
                        contig = convert_chromosome(split_row[0])
                        if contig not in variant_cursors:
                            print('\tCreating variant table for contig {} ...'.
                                  format(contig))
                            new_table = h5file.create_table(
                                where=self.variant_group,
                                name='contig_{}_variants'.format(contig),
                                description=self.VariantTable,
                                title='Variants for contig {}'.format(contig),
                                expectedrows=rows_per_contig[contig])

                            variant_cursors[contig] = new_table.row
                        self._insert_vcf_row(split_row,
                                             variant_cursors[contig])

                for contig in variant_cursors:
                    hdf5 = self.variant_group._f_get_child(
                        'contig_{}_variants'.format(contig))
                    hdf5.flush()
                    # create index
                    for col_name in sorted(
                        ('contig', 'start', 'end', 'pathogenicity',
                         'putative_impact', 'qual', 'rsid')):
                        print('\tGenerating index for column {} in contig {}'.
                              format(col_name, contig))
                        hdf5.cols._f_col(col_name).create_csindex(
                            tmp_dir='/tmp')

        print('\tDone in {}'.format(datetime.datetime.now() - start_time))

    def _insert_vcf_row(self, variant_tuple, variant_cursor):
        """
        Parse a tuple of VCF 4.2 entries and insert them into the current row specified by <variant_cursor>.
        Annotation of CLNSIG field: "Variant Clinical Significance:
            0 - Uncertain significance,
            1 - Not provided,
            2 - Benign,
            3 - Likely benign,
            4 - Likely pathogenic,
            5 - Pathogenic,
            6 - Drug response,
            7 - Histocompatibility, and
            255 - Other".
        """

        # Split up multiple entries in rsid field (if present)
        primary_key_field = list(self.vcf_field_model['primary_key'].keys())[0]
        primary_key_pos, primary_key_parse_func = self.vcf_field_model[
            'primary_key'][primary_key_field]

        primary_keys = primary_key_parse_func(variant_tuple[
            primary_key_pos]).split(';')

        # TODO: also use the 2nd or later RSID
        # if len(primary_keys) > 1:
        #     print('Got multiple ids for same variant: {}'.format(', '.join(primary_keys)))

        for primary_key in primary_keys:
            # Assign the primary key to this variant. Note that re-doing the parsing for each of the ids in a
            # multi-id variant, as we do below, is not the most efficient way of handling this situation.
            # But since this is expected to be rare, we don't optimize yet.
            variant_cursor[primary_key_field] = primary_key

            # Get core fields
            for field_name, (
                    field_position, parse_func
            ) in list(self.vcf_field_model['core_fields'].items()):
                variant_cursor[field_name] = parse_func(variant_tuple[
                    field_position])
            variant_cursor['contig'] = convert_chromosome(variant_cursor[
                'contig'].decode())
            variant_cursor['start'] -= 1

            # end position is with respect to reference sequence
            variant_cursor['end'] = variant_cursor['start'] + len(
                variant_cursor['ref'])
            # get info fields
            split_info = variant_tuple[self.vcf_field_model[
                'info_position']].split(';')

            for field_atom in split_info:
                try:
                    field_key, field_value = field_atom.split('=')
                except ValueError:
                    pass
                else:
                    # Parse the annotation field
                    if field_key == 'ANN':
                        annotations = field_value.split('|')
                        variant_cursor['affected_area'] = annotations[1]
                        variant_cursor['putative_impact'] = annotations[2]
                        variant_cursor['affected_gene'] = annotations[3]

                    if field_key == 'CLNSIG':
                        clinsig_values = [
                            int(v) for v in field_value
                            if field_value in {'2', '3', '4', '5'}
                        ]
                        if len(clinsig_values) > 0:
                            variant_cursor['pathogenicity'] = max(
                                clinsig_values)

                    # Look at the comma-separated _elements of the alt field and the ref field. If any
                    # of them are longer than 1 nucleotide, classify this variant as an indel, otherwise it is a SNP.
                    variant_cursor['is_snp'] = True
                    for allele in variant_cursor['alt'].decode().split(
                            ',') + [variant_cursor['ref']]:
                        if len(allele) > 1:
                            variant_cursor['is_snp'] = False

            # Process "secondary" fields that have field names in column 8 and values in column 9, all
            # colon-delimited.
            try:
                for key, value in zip(
                        variant_tuple[self.vcf_field_model[
                            'secondary_field_keys']].split(':'),
                        variant_tuple[self.vcf_field_model[
                            'secondary_field_values']].split(':')):
                    if key == 'GT':
                        if value == 1 or value == str(1):
                            value = '0/1'
                        variant_cursor[key] = value
            except IndexError:
                print('Error in FORMAT and/or SAMPLE column')

            # update map of rsid to contigs
            self.rsid_index[variant_cursor['rsid'].decode()] = variant_cursor[
                'contig'].decode()

            variant_cursor.append()

    def _tuple_to_dict(variant_tuple, col_names):
        """
        Converts a tuple of variant information as output from a PyTables query, and returns a dictionary
        of field-value pairs. Bytes will be converted to strings as appropriate.
        """

        variant_dict = {}
        for key, value in zip(col_names, variant_tuple):
            try:
                variant_dict[key] = value.decode()
            except AttributeError:
                variant_dict[key] = value
        variant_dict['alleles'] = get_alleles(
            gt_field=variant_dict['GT'],
            ref=variant_dict['ref'],
            alt=variant_dict['alt'])
        return variant_dict

    def get_variant_by_rsid(self, rsid, allow_snps=True, allow_indels=True):
        """
        Return a dictionary of variant attributes for the variant matching <rsid>.
        If <rsid> is not found, return None.
        """

        if rsid in self.rsid_index:
            contig = self.rsid_index[rsid]

            query_string = '(rsid == {})'.format(rsid.encode())

            if allow_snps and allow_indels:
                pass

            elif allow_snps:
                query_string += ' & (is_snp == True)'

            else:
                query_string += ' & (is_snp == False)'

            print('Query: {}'.format(query_string))

            table = self.variant_group._f_get_child(
                'contig_{}_variants'.format(contig))

            query_results = table.read_where(query_string)

            if len(query_results) > 0:
                return self._tuple_to_dict(query_results[0], table.colnames)

        return None

    def get_variants_by_location(self,
                                 contig,
                                 start=0,
                                 end=0,
                                 minimum_quality=0,
                                 allow_snps=True,
                                 allow_indels=False,
                                 minimum_pathogenicity=0,
                                 include_putative_impacts=(),
                                 exclude_putative_impacts=(),
                                 convert_to_dict=False):
        """
        Given a genomic region specified by <contig>, <start>, <end>,
        return all variants overlapping that region.
        :param contig:
        :param start:
        :param end:
        :param minimum_quality:
        :param allow_snps: Whether or not SNPs (variants with either alt or ref fields equal to one nucleotide) should be returned.
        :param allow_indels: Whether or not indels (variants with either alt or ref fields greater than one nucleotide) should be returned.
        :param minimum_pathogenicity: An integer. If set, any variants below this pathogenicity level will be excluded.
        :param include_putative_impacts: An iterable. If set, only variants with these putative impact classifications will be returned
        :param exclude_putative_impacts: An iterable. If set, variants with these putative impact classifications will not be returned
        :param convert_to_dict: If True, return a list of dictionaries of field-value pairs,
         otherwise return a tuple consisting of a list of tuples, and a dictionary mapping field
         names to tuple indices.
        :return:
        """

        start_time = datetime.datetime.now()
        print(
            'Finding all variants of in contig {} ({},{}) with quality > {}, SNPS: {}, Indels: {}...'.
            format(contig, start, end, minimum_quality, allow_snps,
                   allow_indels))
        try:
            table = self.variant_group._f_get_child(
                'contig_{}_variants'.format(contig))
        except tables.exceptions.NoSuchNodeError:
            # We have no variant info for this contig. Either because it's not a valid contig or there were no variants
            # there.
            if convert_to_dict:
                return []
            else:
                return [], {}

        if start > 0 or end > 0:
            query_string = (
                '(contig == {}) & (start >= {}) & (end <= {}) & (qual > {})'.
                format(contig.encode(), start, end, minimum_quality))
        else:
            query_string = ('(contig == {}) & (qual > {})'.format(
                contig.encode(), minimum_quality))

        if allow_snps and allow_indels:
            pass
        elif allow_snps:
            query_string += ' & (is_snp == True)'
        else:
            query_string += ' & (is_snp == False)'

        include_putative_impacts = set(include_putative_impacts)
        exclude_putative_impacts = set(exclude_putative_impacts)
        include_putative_impacts.difference_update(exclude_putative_impacts)

        if include_putative_impacts:
            query_string += ' & ({})'.format(' | '.join([
                '(putative_impact == {})'.format(included_impact.encode())
                for included_impact in include_putative_impacts
            ]))

        if exclude_putative_impacts:
            query_string += ' & ({})'.format(' & '.join([
                '(putative_impact != {})'.format(excluded_impact.encode())
                for excluded_impact in exclude_putative_impacts
            ]))

        if minimum_pathogenicity:
            query_string += ' & (pathogenicity >= {}'.format(
                minimum_pathogenicity)

        print('Query: {}'.format(query_string))

        query_results = table.read_where(query_string)

        print('\tFound {} variants in {}'.format(
            len(query_results), datetime.datetime.now() - start_time))

        if convert_to_dict:
            start_time = datetime.datetime.now()
            print('Converting variants to dictionaries...')
            query_results = [
                self._tuple_to_dict(var, table.colnames)
                for var in query_results
            ]
            print('\tDone in {}'.format(datetime.datetime.now() - start_time))
            return query_results
        else:
            field_mapping = dict([t[::-1] for t in enumerate(table.colnames)])
            return query_results, field_mapping

    def get_variant_columns_by_location(
            self,
            contig,
            start=0,
            end=0,
            minimum_quality=0,
            allow_snps=True,
            allow_indels=False,
            minimum_pathogenicity=0,
            include_putative_impacts=(),
            exclude_putative_impacts=(),
            fields=('start', 'end', 'putative_impact', 'pathogenicity', 'ref',
                    'alt', 'GT')):
        """
        Similar to .get_variants_by_location except that instead of returning an array of structured arrays
        or a list of dictionaries, this method returns a dictionary of arrays, one for each column in the
        result. Iterating over these arrays is much faster than over the row iterator returned by standard
        query methods.
        :param contig:
        :param start:
        :param end:
        :param minimum_quality:
        :param allow_snps:
        :param allow_indels:
        :param minimum_pathogenicity:
        :param include_putative_impacts:
        :param exclude_putative_impacts:
        :param fields: Which fields will be included as columns in the result
        :return:
        """

        start_time = datetime.datetime.now()
        print(
            'Finding all variants of in contig {} ({},{}) with quality > {}, SNPS: {}, Indels: {}...'.
            format(contig, start, end, minimum_quality, allow_snps,
                   allow_indels))

        print('Returning arrays for fields {}.'.format(', '.join(fields)))

        try:
            table = self.variant_group._f_get_child(
                'contig_{}_variants'.format(contig))
        except tables.exceptions.NoSuchNodeError:
            # We have no variant info for this contig. Either because it's not a valid contig or there were no variants
            # there.
            return dict([(field, []) for field in fields])

        if start > 0 or end > 0:
            query_string = (
                '(contig == {}) & (start >= {}) & (end <= {}) & (qual > {})'.
                format(contig.encode(), start, end, minimum_quality))
        else:
            query_string = ('(contig == {}) & (qual > {})'.format(
                contig.encode(), minimum_quality))

        if allow_snps and allow_indels:
            pass
        elif allow_snps:
            query_string += ' & (is_snp == True)'
        else:
            query_string += ' & (is_snp == False)'

        include_putative_impacts = set(include_putative_impacts)
        exclude_putative_impacts = set(exclude_putative_impacts)
        include_putative_impacts.difference_update(exclude_putative_impacts)

        if include_putative_impacts:
            query_string += ' & ({})'.format(' | '.join([
                '(putative_impact == {})'.format(included_impact.encode())
                for included_impact in include_putative_impacts
            ]))

        if exclude_putative_impacts:
            query_string += ' & ({})'.format(' & '.join([
                '(putative_impact != {})'.format(excluded_impact.encode())
                for excluded_impact in exclude_putative_impacts
            ]))

        if minimum_pathogenicity:
            query_string += ' & (pathogenicity >= {}'.format(
                minimum_pathogenicity)

        print('Query: {}'.format(query_string))

        row_coordinates = [r.nrow for r in table.where(query_string)]
        query_results = {}
        for field in fields:
            query_results[field] = table.read_coordinates(
                row_coordinates, field=field)

        print('\tFound {} variants in {}'.format(
            len(query_results[list(query_results.keys())[0]]),
            datetime.datetime.now() - start_time))
        return query_results
