"""Reference path normalization and rotation helpers."""

from hifisr_functions.references import (
    get_rc,
    replace_fasta_id,
    rotate_ref_to_non_repeat_region,
)


FUNCTION_PURITY = {
    "get_rc": "impure",
    "replace_fasta_id": "impure",
    "rotate_ref_to_non_repeat_region": "impure",
}


__all__ = ["get_rc", "replace_fasta_id", "rotate_ref_to_non_repeat_region"]
