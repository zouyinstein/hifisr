# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import polars as pl
import sys
import _bootstrap  # noqa: F401

# Read the excel file
excel_file = sys.argv[1]
engine = "calamine"
out_file = sys.argv[3]
df = pl.read_excel(excel_file, engine=engine)
# df.write_csv(out_file)
df.write_parquet(out_file)
