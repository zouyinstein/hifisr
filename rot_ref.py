import hifisr_functions as hf
import sys


soft_paths_file = sys.argv[1]
soft_paths_dict = { line.split("\t")[0]:line.split("\t")[1] for line in hf.get_file_lines(soft_paths_file) }
genome = sys.argv[2]
old_ref = sys.argv[3]

hf.replace_bait_id(genome, old_ref, genome + "_old.fa")
rot_step = hf.rotate_ref_to_non_repeat_region(genome, genome + "_old.fa", soft_paths_dict, rotation=True) # BLASTn: find large repeats >= 1 kb
