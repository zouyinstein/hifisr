"""Plotting helpers for workflow reports."""

from hifisr_functions.reports import (
    get_gfa_blastn_png,
    get_gfa_reference_images,
    get_gfa_reference_pdf,
    plot_bubble_type_2_rep_raw,
    plot_coverage,
    plot_length_qual,
)


FUNCTION_PURITY = {
    "get_gfa_blastn_png": "impure",
    "get_gfa_reference_images": "impure",
    "get_gfa_reference_pdf": "impure",
    "plot_bubble_type_2_rep_raw": "impure",
    "plot_coverage": "impure",
    "plot_length_qual": "impure",
}


__all__ = [
    "get_gfa_blastn_png",
    "get_gfa_reference_images",
    "get_gfa_reference_pdf",
    "plot_bubble_type_2_rep_raw",
    "plot_coverage",
    "plot_length_qual",
]
