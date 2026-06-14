# hifisr_functions Architecture

This document describes the staged architecture split for `hifisr_functions/`.
The legacy modules remain in place and are still the source of truth for the
current Snakemake workflow and `analysis_scripts/` entry points:

- `base.py`
- `reads.py`
- `references.py`
- `variants.py`
- `annotations.py`
- `transfer.py`
- `reports.py`

The new architecture is currently a compatibility layer. New modules re-export
legacy functions without changing behavior. This allows future code to import
from clearer domain packages while existing scripts continue to work.

Two package names differ from the initial sketch:

- `read_ops/` is used instead of `reads/` so it does not shadow legacy
  `hifisr_functions.reads`.
- `variant_ops/` is used instead of `variants/` so it does not shadow legacy
  `hifisr_functions.variants`.

## Package Map

| New package | Purpose | Legacy sources |
| --- | --- | --- |
| `core/` | Runtime config, command execution, file IO, dependency checks | `base.py` |
| `graph/` | Reference coordinates, sequence paths, GFA/graph visualization, projection to reference | `references.py`, `reports.py` |
| `read_ops/` | Read recruitment, filtering, sampling, correction, stats, coverage reads | `reads.py`, `reports.py`, `variants.py` |
| `variant_ops/` | Junction grouping, SV read alignment, SNV/indel calls, frequency, read linkage | `variants.py`, `annotations.py` |
| `validation/` | Variant flags, NUMT/NUPT transfer-fragment validation, confidence, model selection | `annotations.py`, `transfer.py`, `variants.py`, `references.py`, `reads.py` |
| `report/` | Tables, plots, future HTML reports | `reports.py`, `annotations.py` |
| `workflow/` | Workflow initialization, preflight checkpoints, high-level recipe helpers | `base.py`, `reads.py`, `references.py` |
| `cli/` | Future command-line entry points | package version metadata |

## Import Examples

```python
from hifisr_functions.core.config import load_soft_paths
from hifisr_functions.core.runner import run_checked
from hifisr_functions.graph.coordinate import get_subseq, rotate_fasta
from hifisr_functions.read_ops.recruit import split_mtpt_reads
from hifisr_functions.variant_ops.junction import get_next_groups
from hifisr_functions.validation.numt_nupt import run_transfer_blastn
from hifisr_functions.report.tables import convert_blastn_alignments_to_table
```

Legacy imports remain valid:

```python
import hifisr_functions.reads as hfreads
import hifisr_functions.variants as hfvar
```

## Migration Rule

During migration, add new code to the domain package first, then keep the legacy
module as a wrapper only after all current scripts have moved over. Do not delete
or rename a legacy module while `analysis_scripts/` or `Snakefile` still imports
it directly.
