"""Assembly and polishing helpers used when selecting candidate models."""

from hifisr_functions.reads import correct_reads_by_canu
from hifisr_functions.references import flye_assemble, flye_polish, mecat_cns


FUNCTION_PURITY = {
    "correct_reads_by_canu": "impure",
    "flye_assemble": "impure",
    "flye_polish": "impure",
    "mecat_cns": "impure",
}


__all__ = ["correct_reads_by_canu", "flye_assemble", "flye_polish", "mecat_cns"]
