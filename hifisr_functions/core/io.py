"""File IO helpers for the new architecture."""

from hifisr_functions.base import get_file_lines


FUNCTION_PURITY = {
    "get_file_lines": "impure",
}


__all__ = ["get_file_lines"]
