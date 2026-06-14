"""Shell command runner helpers for the new architecture."""

from hifisr_functions.base import get_cli_output_lines, run_checked


FUNCTION_PURITY = {
    "run_checked": "impure",
    "get_cli_output_lines": "impure",
}


__all__ = ["get_cli_output_lines", "run_checked"]
