import sys

fin = open(sys.argv[1], "rt")
lines = fin.readlines()
for line in lines:
    fields = line.rstrip("\n").split("\t")
    query_id = fields.pop(0)
    length = fields.pop(0)
    type = fields.pop(0)
    count = len(fields)//5
    strand1 = fields.pop(0)
    query_start1 = fields.pop(0)
    query_end1 = fields.pop(0)
    sbjct_start1 = fields.pop(0)
    sbjct_end1 = fields.pop(0)
    strand2 = fields.pop(0)
    query_start2 = fields.pop(0)
    query_end2 = fields.pop(0)
    sbjct_start2 = fields.pop(0)
    sbjct_end2 = fields.pop(0)
    Rep1 = int(query_end1) - int(query_start2) +1
    if Rep1 > 0:
        type = "type_2_rep1"
    elif Rep1 == 0:
        type = "type_2_ref1"
    elif Rep1 < 0:
        type = "type_2_ins1"
    else:
        type = "other"
    new_line = "\t".join((query_id, length, type, str(Rep1), strand1, query_start1, query_end1, sbjct_start1, sbjct_end1, strand2, query_start2, query_end2, sbjct_start2, sbjct_end2))
    print(new_line)

