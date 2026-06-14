"""High-level workflow recipe helpers assembled from legacy routines."""

from hifisr_functions.reads import correct_reads_by_canu
from hifisr_functions.references import (
    aln_to_ref,
    flye_assemble,
    flye_polish,
    mecat_cns,
    simple_draft_asm,
)


FUNCTION_PURITY = {
    "aln_to_ref": "impure",
    "correct_reads_by_canu": "impure",
    "flye_assemble": "impure",
    "flye_polish": "impure",
    "mecat_cns": "impure",
    "simple_draft_asm": "impure",
}


__all__ = [
    "aln_to_ref",
    "correct_reads_by_canu",
    "flye_assemble",
    "flye_polish",
    "mecat_cns",
    "simple_draft_asm",
]
