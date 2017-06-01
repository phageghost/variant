VARIANT_EFFECTS = [
    # Ordered from most to least severe

    # Loss of transcript or exon
    'transcript_ablation',
    'exon_loss_variant',
    # Altered splicing
    'splice_acceptor_variant',
    'splice_donor_variant',
    # Nonsense mutation
    'stop_gained',
    # Frameshift
    'frameshift_variant',
    # Nonstop mutation 1
    'stop_lost',
    # Nonstart mutation
    'start_lost',
    'initiator_codon_variant',
    # Altered transcript 1
    'transcript_amplification',
    'protein_protein_contact',
    'transcript_variant',
    # InDel
    'disruptive_inframe_insertion',
    'disruptive_inframe_deletion',
    'inframe_insertion',
    'inframe_deletion',
    # Altered transcript 2
    'conservative_missense_variant',
    'rare_amino_acid_variant',
    'missense_variant',
    'protein_altering_variant',
    # Altered intragenic region 1
    'splice_region_variant',
    # Nonstop mutation 2
    'incomplete_terminal_codon_variant',
    # Silent mutation
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    # Mutation
    'coding_sequence_variant',
    'exon_variant',
    # Altered miRNA
    'mature_miRNA_variant',
    # Altered 5'UTR
    '5_prime_UTR_variant',
    '5_prime_UTR_premature_start_codon_gain_variant',
    # Altered 3'UTR
    '3_prime_UTR_variant',
    # Altered non-coding exon region
    'non_coding_exon_variant',
    'non_coding_transcript_exon_variant',
    # Altered intragenic region 2
    'intragenic_variant',
    'conserved_intron_variant',
    'intron_variant',
    'INTRAGENIC',
    # Altered nonsense-mediated-decay-target region
    'NMD_transcript_variant',
    # Altered non-coding region
    'non_coding_transcript_variant',
    'nc_transcript_variant',
    # Altered 5'flank site
    'upstream_gene_variant',
    # Altered 3'flank site
    'downstream_gene_variant',
    # Altered transcription-factor-binding region
    'TF_binsing_site_ablation',
    'TFBS_ablation',
    'TF_binding_site_amplification',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'TFBS_variant',
    # Altered regulatory region
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'regulatory_region_variant',
    'regulatory_region',
    'feature_elongation',
    'feature_truncation',
    # Altered intergenic region
    'conserved_intergenic_variant',
    'intergenic_variant',
    'intergenic_region',
    # Others
    'sequence_feature',
]

VCF_COLUMNS = [
    'CHROM',
    'POS',
    'ID',
    'REF',
    'ALT',
    'QUAL',
    'FILTER'
    'INFO',
    'FORMAT',
    # Samples ...
]

CLNSIG_DESCRIPTIONS = {
    0: 'unknown',
    1: 'untested',
    2: 'non-pathogenic',
    3: 'probable-non-pathogenic',
    4: 'probable-pathogenic',
    5: 'pathogenic',
    6: 'drug-response',
    7: 'histocompatibility',
    255: 'other',
}

VCF_ANN_FIELDS = [
    'ALT',
    'effect',
    'impact',
    'gene_name',
    'gene_id',
    'feature_type',
    'feature_id',
    'transcript_biotype',
    'rank',
    'hgvsc',
    'hgvsp',
    'cdna_position',
    'cds_position',
    'protein_position',
    'distance_to_feature',
    'error',
]


# ==============================================================================
#  Function on .VCF information
# ==============================================================================
def parse_vcf_row(vcf_row):
    """
    Parse .VCF row and make a variant dict.
    :param vcf_row: iterable;
    :return: dict; variant dict;
    """

    # CHROM, POS, ID, REF, ALT, QUAL, FILTER
    variant_dict = {
        field: vcf_row[i]
        for (i, field) in enumerate(VCF_COLUMNS[:7])
    }

    # INFO
    without_fields = []  # Some fields are not in field=value format
    for i in vcf_row[7].split(';'):

        if '=' not in i:
            without_fields.append(i)
            continue

        field, value = i.split('=')

        if field == 'ANN':

            # Each ANN is a dict
            ann_dict = {}
            for j, ann in enumerate(value.split(',')):
                ann_split = ann.split('|')
                ann_dict[j] = {
                    VCF_ANN_FIELDS[k]: ann_split[k]
                    for k in range(1, 16)
                }
            variant_dict['ANN'] = ann_dict

        else:
            variant_dict[field] = value

    if without_fields:
        variant_dict['INFO_without_fields'] = '|'.join(without_fields)

    # FORMAT
    format_ = vcf_row[8]
    format_split = format_.split(':')

    # Samples
    if 9 < len(vcf_row):

        # Each sample is a dict
        sample_dict = {}
        for i, sample in enumerate(vcf_row[9:]):
            sample_dict[i] = {
                field: value
                for field, value in zip(format_split, sample.split(':'))
            }
        variant_dict['sample'] = sample_dict

    return variant_dict


def update_vcf_variant_dict(variant_dict):
    """
    Update .VCF variant dict in place.
    :param dict; variant dict
    :return: None
    """

    ref, alt = variant_dict['REF'], variant_dict['ALT']

    variant_dict['variant_type'] = get_variant_type(ref, alt)

    start, end = get_start_and_end_positions(variant_dict['POS'], ref, alt)
    variant_dict['start'] = start
    variant_dict['end'] = end

    variant_dict['confidence'] = 1 - 10**(-float(variant_dict['QUAL']) / 10)

    if 'CLNSIG' in variant_dict:
        variant_dict['clinvar'] = describe_clnsig(variant_dict['CLNSIG'])

    for i, d in variant_dict['ANN'].items():
        d['variant_classification'] = get_variant_classification(d['effect'],
                                                                 ref, alt)

    for i, d in variant_dict['sample'].items():

        if 'GT' in d:
            d['genotype'] = get_vcf_genotype(ref, alt, d['GT'])

        if 'AD' in d and 'DP':
            d['allelic_frequency'] = get_vcf_allelic_frequencies(d['AD'],
                                                                 d['DP'])


def get_vcf_infos(fields, info):
    """
    Get fields from .VCF INFO.
    :param fields: iterable; of str
    :param info: str; .VCF INFO
    :return: list; of str field value
    """

    values = []
    for i in info.split(';'):  # For each INFO

        if '=' not in i:  # Some fields are not in field=value format
            continue

        field, value = i.split('=')

        if field in fields:
            values.append(value)

    return values


def get_vcf_anns(fields, info):
    """
    Get fields from .VCF INFO ANN.
    :param fields: iterable; of str: 'ALT' | 'effect' | 'impact' | 'gene_name'
    | 'gene_id' | 'feature_type' | 'feature_id' | 'transcript_biotype' | 'rank'
    | 'hgvsc' | 'hgvsp' | 'cdna_position' | 'cds_position' | 'protein_position'
    | 'distance_to_feature'| 'error'
    :param info: str; .VCF INFO
    :return: list; of list (ordered by ANN) of str field value
    """

    ann = get_vcf_infos(['ANN'], info).pop()

    anns_values = []

    for a in ann.split(','):
        a_split = a.split('|')

        anns_values.append(
            [a_split[VCF_ANN_FIELDS.index(field)] for field in fields])

    return anns_values


def get_vcf_formats(fields, format_=None, sample=None):
    """
    Get fields' values from .VCF sample.
    :param fields: iterable; of str; .VCF FORMAT fields
    :param format: str; .VCF FORMAT
    :param sample: iterable; of str; .VCF sample
    :return: list; of str (field value)
    """

    format_split = format_.split(':')
    sample_split = sample.split(':')

    return [sample_split[format_split.index(field)] for field in fields]


def get_vcf_genotype(ref, alt, gt=None, format_=None, sample=None):
    """
    Get genotype.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :param gt: str; .VCF sample GT
    :param format_: str; .VCF FORMAT
    :param sample: str; .VCF sample
    :return: str; 'allele_1_sequence|allele_2_sequence'
    """

    gt = gt.replace('/', '|')

    ref_alts = [ref] + alt.split(',')

    return [ref_alts[int(a_gt)] for a_gt in gt.split('|')]


def get_vcf_allelic_frequencies(ad, dp):
    """
    Compute allelic frequency (INFO's AF is a rounded allelic frequency).
    :param ad: str; .VCF sample AD
    :param dp: str; .VCF sample DP
    :return: list; of allelic frequencies
    """

    dp = int(dp)

    return [(int(an_ad) / dp) for an_ad in ad.split(',')]


# ==============================================================================
#  Function on variant information
# ==============================================================================
def get_start_and_end_positions(pos, ref, alt):
    """
    Get variant start and end position.
    :param pos: str; variant position
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: int & int; variant start & end positions
    """

    pos = int(pos)

    if len(ref) == len(alt):
        start, end = pos, pos + len(alt) - 1

    elif len(ref) < len(alt):
        start, end = pos, pos + 1

    else:  # len(alt) < len(ref)
        start, end = pos + 1, pos + len(ref) - len(alt)

    return start, end


def is_inframe(ref, alt):
    """
    Check whether ref-to-alt variant is inframe.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: bool;
    """

    if (len(ref) - len(alt)) % 3:
        return False
    else:
        return True


def get_variant_type(ref, alt):
    """
    Get variant type.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: str; variant type: 'SNP' | 'DNP' | 'TNP' | 'ONP' | 'INS' | 'DEL'
    """

    if len(ref) == len(alt):

        if len(ref) == 1:
            variant_type = 'SNP'

        elif len(ref) == 2:
            variant_type = 'DNP'

        elif len(ref) == 3:
            variant_type = 'TNP'

        else:  # 4 <= len(ref)
            variant_type = 'ONP'

    elif len(ref) < len(alt):
        variant_type = 'INS'

    else:  # len(alt) < len(ref)
        variant_type = 'DEL'

    return variant_type


def describe_clnsig(clnsig, clnsig_descriptions=CLNSIG_DESCRIPTIONS):
    """
    Describe CLNSIG.
    :param clnsig: str; '|' separated: 0 | 1 | 2 | 4 | 5 | 6 | 7 | 255
    :return str; CLNSIG description
    """

    return '|'.join([clnsig_descriptions[int(c)] for c in clnsig.split('|')])


def get_variant_classification(effect, ref, alt):
    """
    Convert .VCF ANN effect to .MAF variant classification.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: str; .MAF variant classification
    """

    variant_type = get_variant_type(ref, alt)

    inframe = is_inframe(ref, alt)

    if effect in (
            'transcript_ablation',
            'exon_loss_variant',
            'splice_acceptor_variant',
            'splice_donor_variant', ):
        variant_classification = 'Splice_Site'

    elif effect in ('stop_gained', ):
        variant_classification = 'Nonsense_Mutation'

    elif variant_type == 'INS' and (effect == 'frameshift_variant' or
                                    (not inframe and effect in (
                                        'protein_protein_contact',
                                        'protein_altering_variant', ))):
        variant_classification = 'Frame_Shift_Ins'

    elif variant_type == 'DEL' and (effect == 'frameshift_variant' or
                                    (not inframe and effect in (
                                        'protein_protein_contact',
                                        'protein_altering_variant', ))):
        variant_classification = 'Frame_Shift_Del'

    elif effect in ('stop_lost', ):
        variant_classification = 'Nonstop_Mutation'

    elif effect in (
            'start_lost',
            'initiator_codon_variant', ):
        variant_classification = 'Translation_Start_Site'

    elif variant_type == 'INS' and inframe and effect in (
            'protein_protein_contact',
            'disruptive_inframe_insertion',
            'inframe_insertion',
            'protein_altering_variant', ):
        variant_classification = 'In_Frame_Ins'

    elif variant_type == 'DEL' and inframe and effect in (
            'protein_protein_contact',
            'disruptive_inframe_deletion',
            'inframe_deletion',
            'protein_altering_variant', ):
        variant_classification = 'In_Frame_Del'

    elif effect in (
            'transcript_variant',
            'conservative_missense_variant',
            'rare_amino_acid_variant',
            'missense_variant',
            'coding_sequence_variant', ) or (
                variant_type not in ('INS', 'DEL') and
                effect == 'protein_protein_contact'):
        variant_classification = 'Missense_Mutation'

    elif effect in (
            'transcript_amplification',
            'splice_region_variant',
            'intragenic_variant',
            'conserved_intron_variant',
            'intron_variant',
            'INTRAGENIC', ):
        variant_classification = 'Intron'

    elif effect in (
            'incomplete_terminal_codon_variant',
            'start_retained_variant',
            'stop_retained_variant',
            'synonymous_variant',
            'NMD_transcript_variant', ):
        variant_classification = 'Silent'

    elif effect in (
            'exon_variant',
            'mature_miRNA_variant',
            'non_coding_exon_variant',
            'non_coding_transcript_exon_variant',
            'non_coding_transcript_variant',
            'nc_transcript_variant', ):
        variant_classification = 'RNA'

    elif effect in (
            '5_prime_UTR_variant',
            '5_prime_UTR_premature_start_codon_gain_variant', ):
        variant_classification = '5\'UTR'

    elif effect in ('3_prime_UTR_variant', ):
        variant_classification = '3\'UTR'

    elif effect in (
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
        variant_classification = 'IGR'

    elif effect in ('upstream_gene_variant', ):
        variant_classification = '5\'Flank'

    elif effect in ('downstream_gene_variant', ):
        variant_classification = '3\'Flank'

    elif effect in ('sequence_feature', ):
        variant_classification = 'Targeted_Region'

    else:
        print('Unknown: effect={} & variant_type={} & inframe={}.'.format(
            effect, variant_type, inframe))
        variant_classification = 'Targeted_Region'

    return variant_classification
