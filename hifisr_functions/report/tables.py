"""Tabular report and annotation-output helpers."""

from hifisr_functions.annotations import (
    add_depth_and_frq,
    combine_variant_anno,
    get_variant_types,
)
from hifisr_functions.reports import convert_blastn_alignments_to_table


FUNCTION_PURITY = {
    "add_depth_and_frq": "impure",
    "combine_variant_anno": "impure",
    "convert_blastn_alignments_to_table": "impure",
    "get_variant_types": "impure",
}


__all__ = [
    "add_depth_and_frq",
    "combine_variant_anno",
    "convert_blastn_alignments_to_table",
    "get_variant_types",
]
