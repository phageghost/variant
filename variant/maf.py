from numpy import argmin, empty
from pandas import DataFrame, read_csv

from .variant import MUTATION_EFFECT_RANKING


def convert_vcf_to_maf(vcf,
                       sample_name,
                       ensg_to_entrez,
                       n_anns=1,
                       maf_file_path=None):
    """
    .VCF ==> .MAF.
    :param vcf: DataFrame; .VCF data
    :param sample_name:
    :param ensg_to_entrez: str; File path to a id mapping file (ENSG<\t>Entrez)
    :param n_anns: int;
    :param maf_file_path: str; .MAF file path.
    :return: DataFrame; .MAF data
    """

    # Read ENSG-to-Entrez dict
    ensg_to_entrez_dict = read_csv(ensg_to_entrez, index_col=0).to_dict()

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
        if i % 1000 == 0:
            print('\t@ {} ...'.format(i + 1))

        chrom, pos, id_, ref, alt, info = vcf_row.iloc[[0, 1, 2, 3, 4, 7]]

        start, end = _get_variant_start_and_end_positions(pos, ref, alt)

        if (len(ref) - len(alt)) % 3:
            inframe = False
        else:
            inframe = True

        vt = get_variant_type(ref, alt)

        for i_s in info.split(';'):  # For each INFO
            if i_s.startswith('ANN='):  # ANN
                # Strip 'ANN=' prefix
                i_s = i_s[4:]
                anns = i_s.split(',')

                for a in anns[:n_anns]:  # For each ANN
                    a_s = a.split('|')

                    effect = a_s[1]

                    vc = get_maf_variant_classification(effect, vt, inframe)

                    gene_name = a_s[3]

                    gene_id = ensg_to_entrez_dict.get(a_s[4], 0)

        tmp[i] = gene_name, gene_id, chrom, start, end, vc, vt, id_, ref, alt

    maf[[
        'Hugo_Symbol', 'Entrez_Gene_Id', 'Chromosome', 'Start_Position',
        'End_Position', 'Variant_Classification', 'Variant_Type', 'dbSNP_RS',
        'Reference_Allele', 'Tumor_Seq_Allele1'
    ]] = tmp
    maf['Strand'] = '+'
    maf[[
        'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
        'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID'
    ]] = sample_name

    # Save
    if maf_file_path:
        if not maf_file_path.endswith('.maf'):
            maf_file_path += '.maf'
        maf.to_csv(maf_file_path, sep='\t', index=None)

    return maf


def get_maf_variant_classification(es, vt, inframe):
    """

    :param es: str; effect or effects concatenated by '&'
    :param vt: str; Variant type
    :param inframe: bool;
    :return: str; MAF variant classification
    """

    # Split effects
    es = es.split('&')

    # Return the worst effect
    return convert_ann_effect_to_maf_variant_classification(
        es[argmin([MUTATION_EFFECT_RANKING.index(e) for e in es])], vt,
        inframe)


def convert_ann_effect_to_maf_variant_classification(e, vt, inframe):
    """
    Convert .VCF's variant classification to .MAF's.
    :param e: str; .VCF effect
    :param vt: str; variant type; 'SNP' | 'DNP' | 'TNP' | 'ONP' | 'INS' | 'DEL'
    :param inframe: bool; inframe event or not
    :return: str; .MAF variant classification
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
    """
    Convert .MAF's variant classification to MUTSIG's.
    :param vc: str; MAF variant classification
    :return: str; 'noncoding' | 'null' | 'silent' | 'nonsilent' | 'effect'
    """

    return {
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
    }[vc]
