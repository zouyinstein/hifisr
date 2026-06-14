"""Confidence-support helpers based on depth and read-frequency evidence."""

from hifisr_functions.annotations import add_depth_and_frq
from hifisr_functions.variants import check_FL_and_multi


FUNCTION_PURITY = {
    "add_depth_and_frq": "impure",
    "check_FL_and_multi": "impure",
}


__all__ = ["add_depth_and_frq", "check_FL_and_multi"]
