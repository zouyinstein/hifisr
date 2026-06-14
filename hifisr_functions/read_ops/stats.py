"""Read statistics, quality plotting, and error-rate helpers."""

from hifisr_functions.reads import calc_error_rate
from hifisr_functions.reports import get_fastq_stats, plot_length_qual


FUNCTION_PURITY = {
    "calc_error_rate": "impure",
    "get_fastq_stats": "impure",
    "plot_length_qual": "impure",
}


__all__ = ["calc_error_rate", "get_fastq_stats", "plot_length_qual"]
