# HiFiSR minimap2 structural-read alignment note

Date: 2026-06-14

## Scope

This note records the current minimap2 experiment for the W3-5-2 compatibility
run. The legacy BLASTN path remains the default.

Baseline:

- `/Users/zouyinstein-m4max/Documents/Codex/hifisr-dev/hifisr/results/W3-5-2`

minimap2 output:

- `/Users/zouyinstein-m4max/Documents/Codex/hifisr-dev/hifisr/results/W3-5-2_minimap2`

Input reads:

- `/Users/zouyinstein-m4max/Documents/Codex/hifisr-dev/data/W3-5-2.fastq.gz`

## Implementation

The legacy BLASTN functions below are intentionally unchanged:

- `hifisr_functions.variants.run_blastn_sorter_single`
- `hifisr_functions.variants.run_multi_threads_blastn`
- `hifisr_functions/references.py`

The minimap2 path now lives in `hifisr_functions/variants.py` under distinct
function names:

- `run_minimap2_sorter_single`
- `run_multi_threads_minimap2`

`analysis_scripts/get_variants_in_reads.py` keeps the original BLASTN call as a
comment and currently calls `hfvar.run_multi_threads_minimap2(...)` to generate
the same legacy `all_sorted_blastn_alignments.txt` filename and internal format.
The downstream parsing, grouping, plotting, coverage, and bcftools calls remain
unchanged.

`references.py` still uses BLASTN, including reference repeat/rotation and
`aln_to_ref`. This keeps repeat identification separate from the read-level SV
alignment replacement until a minimap2 replacement can be evaluated with a
truth-independent rule.

## Current minimap2 command

The active read-level structural alignment parameters are:

```text
minimap2 -t 1 -x map-hifi -c -k 11 -w 7 --secondary=yes -N 200 -p 0.01 <ref.fa> <read.fa>
```

The minimap2 PAF output is converted into the existing BLASTN-compatible
`all_sorted_blastn_alignments.txt` internal format. File names and downstream
parsers remain unchanged.

Before writing the BLASTN-compatible summary, the minimap2 path applies a
deterministic terminal micro-indel normalization step. This step uses only the
read sequence, reference sequence, and minimap2 PAF coordinates:

- inspect at most 8 bp outside each alignment end
- consider only alignments at least 500 bp long
- allow a terminal extension only when one side has a 1-2 bp indel and the
  remaining bases match exactly
- choose the extension with the most matched bases, and for ties choose the
  shortest endpoint shift

This is not a BLASTN-guided selection rule; it can run on new samples without an
answer key. The sidecar `whole_read_evidence.tsv` remains raw minimap2 evidence.

## Sidecar Evidence

Each minimap2 run also writes:

- `backup_info/whole_read_evidence.tsv`
- `backup_info/minimap2_runtime.tsv`
- `backup_info/minimap2_runtime_summary.tsv`

Columns:

- `read_id`
- `read_length`
- `query_start`
- `query_end`
- `target_id`
- `target_start`
- `target_end`
- `strand`
- `identity`
- `mapq`
- `alignment_role`
- `cigar`
- `crosses_junction`
- `junction_id`
- `crosses_repeat_choice`
- `repeat_choice_id`
- `source_compat_file`

The junction/repeat-choice columns are present as sidecar placeholders and are
currently `not_evaluated`; they should be filled by a later graph-aware evidence
projection step rather than by the BLASTN-compatibility summary.

`minimap2_runtime.tsv` records per-read minimap2 and terminal-extension timing.
`minimap2_runtime_summary.tsv` records the minimap2 command options, terminal
extension parameters, worker wall time for generating
`all_sorted_blastn_alignments.txt`, summed minimap2 time, and summed terminal
extension time.

## W3-5-2 Full-Run Result

Final comparison report:

- `/Users/zouyinstein-m4max/Documents/Codex/hifisr-dev/hifisr/results/W3-5-2_minimap2/comparison/minimap2_vs_blastn_terminal_microindel_report.md`

Summary against the BLASTN baseline:

| Run | shared reads | exact line match | `aln_type` match | mean abs percent diff | reads with >1 percent diff |
| --- | ---: | ---: | ---: | ---: | ---: |
| mito/run_2 | 4000 | 2773 / 4000 | 3750 / 4000 (93.75%) | 0.1528 | 60 |
| plastid/run_2 | 4000 | 2775 / 4000 | 3972 / 4000 (99.30%) | 0.0113 | 5 |
| mito/run_3 | 4000 | 2773 / 4000 | 3750 / 4000 (93.75%) | 0.1528 | 60 |

Filtered variant-table row counts match the baseline:

- mito/run_2: old 3, new 3
- plastid/run_2: old 0, new 0
- mito/run_3: old 1, new 1

For mito, `total_count`, `depth`, and `frequency` match on the shared variant
rows, but multi-allelic `alt` string order can differ. The remaining table
differences are therefore mostly allele-order representation rather than count
or depth changes.

## Parameter Notes

The first full run used `-x map-hifi -c --secondary=yes -N 50 -p 0.5`. The
current parameters improved mito `aln_type` matching from about 89.7% to 92.65%
and reduced reads with >1 percent-total difference from 148 to 60.

A more aggressive probe, `-k 9 -w 5 -N 300 -p 0.001 -m 20 -n 2 -s 40`, improved
the mito probe further but produced many more PAF records and was too slow on
the plastid probe to use as the current default. It remains a candidate for a
future explicit high-sensitivity mode.

Additional boundary-scoring probes showed that low gap-open settings can make
some high-count `type_2_subtype_rep_NA` coordinates and `mid_olp_1` values match
the BLASTN baseline more closely on W3-5-2. Those probes are diagnostic only:
the production minimap2 path must not select between multiple minimap2 profiles
by retrospective closeness to BLASTN, because a real new sample will not have a
BLASTN answer key.

The accepted optimization is the single-minimap2 terminal micro-indel
normalization above. On the mito/run_2 lightweight summary probe, using
`window=8`, `max_gap=2`, and shortest-shift tie-breaking brought the high
`subgroup_count` `type_2_subtype_rep_NA` groups back to the same `se1/ss2`
coordinates and nearly the same `mid_olp_1` values as the BLASTN baseline. The
remaining large differences are the low `olp=4` groups, where reads are
classified into different alignment types rather than merely shifted by a few
endpoint bases.

Manual biological interpretation for W3-5-2: the BLASTN `olp=42/43` near-diagonal
signal is a 29 bp tandem-repeat copy deletion caused by nuclear NUMT
interference. The BLASTN `olp=4` near-diagonal groups similarly reflect short
nearby deletions, such as 48 bp or 70 bp, from NUMT interference. Moving these
reads into `type_1` in the minimap2 path is acceptable because the downstream
bcftools step can mine the deletion evidence without treating it as a structural
repeat group.

## Next Recommendation

Add an explicit `SV_aligner` or `SV_aligner_mode` parameter later, defaulting to
`blastn`. Suggested values:

- `blastn`: current legacy behavior
- `minimap2`: current compatibility mode
- `minimap2_sensitive`: optional high-sensitivity mode after runtime testing
