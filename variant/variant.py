from re import split

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
def parse_vcf_row(vcf_row, n_anns=1):
    """
    Parse .VCF row and make a variant dict.
    :param vcf_row: iterable;
    :param n_anns: int; number of ANNs to parse
    :return: dict; variant dict; won't have missing non-main field
    """

    # TODO: Elegantify
    variant_dict = {
        'CHROM': cast_vcf_field_value('CHROM', vcf_row[0]),
        'POS': cast_vcf_field_value('POS', vcf_row[1]),
        'ID': cast_vcf_field_value('ID', vcf_row[2]),
        'REF': cast_vcf_field_value('REF', vcf_row[3]),
        'ALT': cast_vcf_field_value('ALT', vcf_row[4]),
        'QUAL': cast_vcf_field_value('QUAL', vcf_row[5]),
        'FILTER': cast_vcf_field_value('FILTER', vcf_row[6]),
    }

    ref = variant_dict['REF']
    alt = variant_dict['ALT']

    # Variant type
    if alt and alt != '.':
        variant_dict['variant_type'] = get_variant_type(variant_dict['REF'],
                                                        alt)

    # Samples
    variant_dict['samples'] = []
    format_split = vcf_row[8].split(':')
    for i, s in enumerate(vcf_row[9:]):
        s_d = {'sample_id': i + 1}

        # Sample
        for format_field, sample_value in zip(format_split, s.split(':')):
            s_d[format_field] = cast_vcf_field_value(format_field,
                                                     sample_value)

        # Genotype
        if alt and alt != '.':
            ref_alts = [ref] + alt.split(',')
            s_d['genotype'] = [ref_alts[int(gt)] for gt in s_d['GT']]
        else:
            s_d['genotype'] = [ref] * 2

        # Allelic frequency
        s_d['allelic_frequency'] = get_allelic_frequencies(vcf_row[8], s)

        variant_dict['samples'].append(s_d)

    info_split = vcf_row[7].split(';')
    for i_s in info_split:
        if i_s.startswith('ANN='):
            anns = {}
            for i, a in enumerate(i_s.split(',')[:n_anns]):
                a_split = a.split('|')
                anns[i] = {ANN_FIELDS[j]: a_split[j] for j in range(1, 15)}
            variant_dict['ANN'] = anns
        else:
            field, value = i_s.split('=')
            variant_dict[field] = cast_vcf_field_value(field, value)

    return variant_dict


def get_info(vcf_row, field):
    """
    Get field from variant INFO.
    :param vcf_row: iterable; a .VCF row
    :param field: str;
    :return: str; field value
    """

    for fv in vcf_row[7].split(';'):  # For each INFO
        f, v = fv.split('=')
        if f == field:
            return v


def get_ann(vcf_row, field, n_anns=1, ann_fields=ANN_FIELDS):
    """
    Get field from variant ANN(s).
    :param vcf_row: iterable; a .VCF row
    :param field: str; 'ALT' | 'effect' | 'impact' | 'gene_name' | 'gene_id' |
    'feature_type' | 'feature_id' | 'transcript_biotype' | 'rank' | 'hgvsc' |
    'hgvsp' | 'cdna_position' | 'cds_position' | 'protein_position' |
    'distance_to_feature'| 'error'
    :param n_anns: int; number of ANNs to get field from
    :return: list; of field values
    """

    # ANN is in INFO, which is the 7th .VCF column
    for f in vcf_row[7].split(';'):  # For each INFO

        if f.startswith('ANN='):  # ANN

            # Strip 'ANN=' prefix
            f = f[4:]

            vs = []

            # 1 variant can have multiple ANNs split by ','
            anns = f.split(',')

            for a in anns[:n_anns]:  # For each ANN
                vs.append(a.split('|')[ann_fields.index(field)])

            return vs


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
            'ID': str,
            'QUAL': float,
            'GT': lambda x: split('[|/]', x),
            'AD': lambda x: x.split(','),
            'CLNSIG': lambda x: max([int(s) for s in split('[,|]', x)]),
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
