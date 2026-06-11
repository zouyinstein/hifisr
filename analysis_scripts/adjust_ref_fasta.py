# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

import sys
import _bootstrap  # noqa: F401
import hifisr_functions.base as hfbase
import hifisr_functions.references as hfref


# Usage: python adjust_ref_fasta.py soft_paths_file genome input_ref_fasta
# Load the soft paths
soft_paths_file = sys.argv[1]
soft_paths_dict = hfbase.load_soft_paths(soft_paths_file)

# Parse other arguments
genome = sys.argv[2] # mito or plastid
input_ref_fasta = sys.argv[3]
hfref.replace_fasta_id(genome, input_ref_fasta, "input_" + genome + ".fa")

# Rotate to the middle of the longest non-repeat region
rot_step = hfref.rotate_ref_to_non_repeat_region(genome, "input_" + genome + ".fa", soft_paths_dict, rotation=True)
# or rotate by a fixed step clockwise
# rot_step = 1000
# hfref.rotate_fasta("input_" + genome + ".fa", genome + "_rotated_" + str(rot_step) + ".fasta", rot_step)
