"""Read recruitment and reference-binning helpers."""

from hifisr_functions.reads import (
    cal_ID_coverage,
    replace_reads_id,
    split_mtpt_reads,
    split_reads_by_contig,
)
from hifisr_functions.variants import get_cov_reads


FUNCTION_PURITY = {
    "cal_ID_coverage": "impure",
    "get_cov_reads": "impure",
    "replace_reads_id": "impure",
    "split_mtpt_reads": "impure",
    "split_reads_by_contig": "impure",
}


__all__ = [
    "cal_ID_coverage",
    "get_cov_reads",
    "replace_reads_id",
    "split_mtpt_reads",
    "split_reads_by_contig",
]
