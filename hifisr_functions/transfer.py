# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

"""Transfer-fragment analysis helpers.

This module provides a stable `hftrans` import alias for transfer-fragment
workflows. The implementation currently delegates to routines maintained in
`hifisr_functions.annotations`.
"""

from hifisr_functions.annotations import (
    get_raw_transfer_groups,
    merge_fragments_files,
    merge_numt_nupt,
    remove_inner_fragments,
    run_transfer_blastn,
)


# Function purity marker:
# pure = deterministic from explicit inputs with no file, shell, environment,
# logging, or input-mutation side effects; impure = otherwise.
FUNCTION_PURITY = {
    "get_raw_transfer_groups": "impure",
    "merge_fragments_files": "impure",
    "merge_numt_nupt": "impure",
    "remove_inner_fragments": "impure",
    "run_transfer_blastn": "impure",
}


__all__ = [
    "get_raw_transfer_groups",
    "merge_fragments_files",
    "merge_numt_nupt",
    "remove_inner_fragments",
    "run_transfer_blastn",
]
