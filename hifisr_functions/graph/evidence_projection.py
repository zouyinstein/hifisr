"""Helpers that project polished assembly evidence back to references."""

from hifisr_functions.references import aln_to_ref


FUNCTION_PURITY = {
    "aln_to_ref": "impure",
}


__all__ = ["aln_to_ref"]
