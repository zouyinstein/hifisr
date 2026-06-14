"""Variant flagging and annotation helpers."""

from hifisr_functions.annotations import combine_variant_anno, get_variant_types
from hifisr_functions.variants import snv_or_indel


FUNCTION_PURITY = {
    "combine_variant_anno": "impure",
    "get_variant_types": "impure",
    "snv_or_indel": "impure",
}


__all__ = ["combine_variant_anno", "get_variant_types", "snv_or_indel"]
