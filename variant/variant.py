from re import split

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
]

ANN_FIELDS = [
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

MUTATION_EFFECT_RANKING = [
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


# ==============================================================================
#  Function on .VCF row
# ==============================================================================
def parse_vcf_row(vcf_row):
    """
    Parse .VCF row and make a variant dict.
    :param vcf_row: iterable;
    :return: dict; variant dict; won't have missing non-main field
    """

    # CHROM, ID, REF, ALT, QUAL, FILTER
    variant_dict = {
        field: cast_vcf_field_value(field, vcf_row[i])
        for (i, field) in enumerate(VCF_COLUMNS[:7])
    }

    ref = variant_dict['REF']
    alt = variant_dict['ALT']

    # Variant type
    variant_dict['variant_type'] = get_variant_type(ref, alt)

    # INFO
    for info in vcf_row[7].split(';'):
        field, value = info.split('=')

        if field == 'ANN':
            # Use 1st ANN
            a = value.split(',')[0].split('|')
            variant_dict['ANN'] = {ANN_FIELDS[i]: a[i] for i in range(1, 15)}

        else:
            variant_dict[field] = cast_vcf_field_value(field, value)

    # FORMAT
    format_ = vcf_row[8]
    format_split = format_.split(':')

    # Samples
    variant_dict['samples'] = []
    for i, sample in enumerate(vcf_row[9:]):

        # Each sample is a dict
        sample_dict = {'sample_id': i + 1}

        for field, value in zip(format_split, sample.split(':')):
            sample_dict[field] = cast_vcf_field_value(field, value)

        # Genotype
        ref_alts = [ref] + alt.split(',')
        sample_dict[
            'genotype'] = [ref_alts[int(gt)] for gt in sample_dict['GT']]

        # Allelic frequency
        sample_dict['allelic_frequency'] = get_allelic_frequencies(format_,
                                                                   sample)

        variant_dict['samples'].append(sample_dict)

    return variant_dict


def get_info(field, vcf_row=None, info=None):
    """
    Get field from variant INFO.
    :param field: str;
    :param vcf_row: iterable; a .VCF row
    :param info: str; INFO
    :return: str; field value
    """

    if not info:
        info = vcf_row[7]

    for fv in info.split(';'):  # For each INFO

        f, v = fv.split('=')

        if f == field:
            return v


def get_ann(fields, vcf_row=None, info=None):
    """
    Get field from variant ANN.
    :param fields: iterable; of str: 'ALT' | 'effect' | 'impact' | 'gene_name'
    | 'gene_id' | 'feature_type' | 'feature_id' | 'transcript_biotype' | 'rank'
    | 'hgvsc' | 'hgvsp' | 'cdna_position' | 'cds_position' | 'protein_position'
    | 'distance_to_feature'| 'error'
    :param vcf_row: iterable; a .VCF row
    :param info: str; INFO
    :return: list; of str field value
    """

    # ANN is in INFO, which is the 7th .VCF column
    if not info:
        info = vcf_row[7]

    for fv in info.split(';'):  # For each INFO

        if '=' not in fv:  # Some fields are not in field=value format
            # print('{} not in FIELD=VALUE in INFO.'.format(fv))
            continue

        f, v = fv.split('=')
        if f == 'ANN':

            # Variant can have multiple ANNs, but use the 1st ANN
            a_split = v.split(',')[0].split('|')

            return [a_split[ANN_FIELDS.index(f)] for f in fields]


# ==============================================================================
#  Function on .VCF row element
# ==============================================================================
def cast_vcf_field_value(field, value):
    """
    Cast .VCF field's value.
    :param field: str; field
    :param value: str; value
    :return: int | float | tuple;
    """

    try:
        return {
            'POS': int,
            'QUAL': float,
            'GT': lambda v: split('[|/]', v),
            'AD': lambda v: v.split(','),
            'CLNSIG': lambda v: max([int(s) for s in split('[,|]', v)]),
        }[field](value)

    except ValueError:
        return value


def get_variant_start_and_end_positions(pos, ref, alt):
    """
    Get variant start and end position.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: tuple; of ints; (start, end)
    """

    if len(ref) == len(alt):
        s, e = pos, pos + len(alt) - 1

    elif len(ref) < len(alt):
        s, e = pos, pos + 1

    else:  # len(alt) < len(ref)
        s, e = pos + 1, pos + len(ref) - len(alt)

    return s, e


def get_variant_type(ref, alt):
    """
    Get variant type.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: str; variant type; 'SNP' | 'DNP' | 'TNP' | 'ONP' | 'INS' | 'DEL'
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


def get_genotype(format_, sample):
    """
    Get genotype.
    :param format_: str; .VCF FORMAT column
    :param sample: str; .VCF sample column
    :return: str;
    """
    return


def get_allelic_frequencies(format_, sample):
    """
    Compute allelic frequency (INFO's AF is a rounded allelic frequency).
    :param format_: str; .VCF FORMAT column
    :param sample: str; .VCF sample column
    :return: list; of allelic frequencies
    """

    format_split = format_.split(':')
    sample_split = sample.split(':')

    ad_i = format_split.index('AD')
    dp_i = format_split.index('DP')

    ads = [int(ad) for ad in sample_split[ad_i]]
    dp = int(sample_split[dp_i])

    return ['{:.3f}'.format(ad / dp) for ad in ads]
