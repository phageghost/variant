from numpy import empty
from pandas import DataFrame, read_csv

from .variant import (get_start_and_end_positions, get_variant_classification,
                      get_variant_type, get_vcf_anns, is_inframe)


def make_maf_from_vcf(vcf,
                      ensg_to_entrez,
                      sample_name='Sample',
                      maf_file_path=None):
    """
    Make .MAF from .VCF.
    :param vcf: DataFrame; .VCF data
    :param ensg_to_entrez: str; File path to a ID-mapping file (ENSG<\t>Entrez)
    :param sample_name: str;
    :param maf_file_path: str; .MAF file path.
    :return: DataFrame; .MAF data
    """

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

    # Read ENSG-to-Entrez dict
    ensg_to_entrez_dict = read_csv(ensg_to_entrez, index_col=0).to_dict()

    print('Iterating through VCF rows ...')
    tmp = empty((vcf.shape[0], 10), dtype=object)
    for i, vcf_row in vcf.iterrows():  # For each VCF row
        if i % 1000 == 0:
            print('\t@ {} ...'.format(i + 1))

        chrom, pos, id_, ref, alt, info = vcf_row[[0, 1, 2, 3, 4, 7]]

        start, end = get_start_and_end_positions(pos, ref, alt)

        inframe = is_inframe(ref, alt)

        variant_type = get_variant_type(ref, alt)

        effect, gene_name, gene_id = get_vcf_anns(
            ['effect', 'gene_name', 'gene_id'], info=info)

        entrez_gene_id = ensg_to_entrez_dict.get(gene_id)

        variant_classification = get_variant_classification(
            effect, variant_type, inframe)

        tmp[i] = gene_name, entrez_gene_id, chrom, start, end, variant_classification, variant_type, id_, ref, alt

    maf[[
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
    maf['Strand'] = '+'
    maf[[
        'Tumor_Sample_Barcode',
        'Matched_Norm_Sample_Barcode',
        'Tumor_Sample_UUID',
        'Matched_Norm_Sample_UUID',
    ]] = sample_name

    # Save
    if maf_file_path:
        if not maf_file_path.endswith('.maf'):
            maf_file_path += '.maf'
        maf.to_csv(maf_file_path, sep='\t', index=None)

    return maf


def get_mutsig_effect(variant_classification):
    """
    Convert .MAF variant classification to MUTSIG effect.
    :param variant_classification: str; .MAF variant classification
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
    }[variant_classification]
