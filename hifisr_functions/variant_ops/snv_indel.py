"""SNV and indel calling helpers."""

from hifisr_functions.variants import (
    run_bcftools,
    run_multi_threads_bcftools,
    snv_or_indel,
)


FUNCTION_PURITY = {
    "run_bcftools": "impure",
    "run_multi_threads_bcftools": "impure",
    "snv_or_indel": "impure",
}


__all__ = ["run_bcftools", "run_multi_threads_bcftools", "snv_or_indel"]
