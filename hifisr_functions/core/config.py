"""Configuration and runtime dependency helpers.

This module is the new-architecture entry point for configuration-related
helpers. The current implementation delegates to the legacy `base.py` module so
existing workflow behavior remains unchanged.
"""

from hifisr_functions.base import (
    PYTHON_IMPORT_NAMES,
    REQUIRED_WORKFLOW_TOOL_GROUPS,
    REQUIRED_WORKFLOW_TOOLS,
    check_python_packages,
    load_soft_paths,
    parse_python_requirements,
    require_soft_paths,
    validate_soft_paths,
)


FUNCTION_PURITY = {
    "load_soft_paths": "impure",
    "validate_soft_paths": "impure",
    "require_soft_paths": "impure",
    "parse_python_requirements": "impure",
    "check_python_packages": "impure",
}


__all__ = [
    "PYTHON_IMPORT_NAMES",
    "REQUIRED_WORKFLOW_TOOL_GROUPS",
    "REQUIRED_WORKFLOW_TOOLS",
    "check_python_packages",
    "load_soft_paths",
    "parse_python_requirements",
    "require_soft_paths",
    "validate_soft_paths",
]
