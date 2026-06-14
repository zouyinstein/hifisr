"""Read filtering, sampling, and correction helpers."""

from hifisr_functions.reads import (
    correct_reads_by_canu,
    filt_length_qual,
    random_sampling,
)


FUNCTION_PURITY = {
    "correct_reads_by_canu": "impure",
    "filt_length_qual": "impure",
    "random_sampling": "impure",
}


__all__ = ["correct_reads_by_canu", "filt_length_qual", "random_sampling"]
