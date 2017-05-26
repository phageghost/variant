from tabix import open

from .variant import parse_vcf_row


def get_variants_by_tabix(sample_vcf,
                          chrom=None,
                          start=None,
                          end=None,
                          query_str=None,
                          reference_vcf=None):
    """
    Get .VCF variants by tabix.
    :param sample_vcf: str or pytabix handler;
    :param chrom: str;
    :param start: int;
    :param end: int;
    :param query_str: str; 'chr:start-end'
    :param reference_vcf: str or pytabix handler;
    :return: list; of variant dicts
    """

    if isinstance(sample_vcf, str):  # Open sample VCF
        sample_vcf = open(sample_vcf)

    # Query sample
    if query_str:
        records = sample_vcf.querys(query_str)
    else:
        records = sample_vcf.query(chrom, start - 1, end)

    if reference_vcf and len(
            list(records)
    ) == 0:  # If reference VCF is available and querying sample failed

        if isinstance(reference_vcf, str):  # Open reference VCF
            reference_vcf = open(reference_vcf)

        # Query reference
        if query_str:
            records = reference_vcf.querys(query_str)
        else:
            records = reference_vcf.query(chrom, start - 1, end)

    return [parse_vcf_row(r) for r in records]
