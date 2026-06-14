"""Workflow initialization helpers."""

from hifisr_functions.base import load_soft_paths, parse_python_requirements


FUNCTION_PURITY = {
    "load_soft_paths": "impure",
    "parse_python_requirements": "impure",
}


__all__ = ["load_soft_paths", "parse_python_requirements"]
