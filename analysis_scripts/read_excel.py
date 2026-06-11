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
