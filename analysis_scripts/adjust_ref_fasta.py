import sys
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
