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
    strand3 = fields.pop(0)
    query_start3 = fields.pop(0)
    query_end3 = fields.pop(0)
    sbjct_start3 = fields.pop(0)
    sbjct_end3 = fields.pop(0)
    strand4 = fields.pop(0)
    query_start4 = fields.pop(0)
    query_end4 = fields.pop(0)
    sbjct_start4 = fields.pop(0)
    sbjct_end4 = fields.pop(0)
    strand5 = fields.pop(0)
    query_start5 = fields.pop(0)
    query_end5 = fields.pop(0)
    sbjct_start5 = fields.pop(0)
    sbjct_end5 = fields.pop(0)

    Rep1 = int(query_end1) - int(query_start2) +1
    Rep2 = int(query_end2) - int(query_start3) +1
    Rep3 = int(query_end3) - int(query_start4) +1
    Rep4 = int(query_end4) - int(query_start5) +1
    if Rep1 > 0 and Rep2 > 0:
        type = "type_5_rep1_rep2"
    elif Rep1 > 0 and Rep2 == 0:
        type = "type_5_rep1_ref2"
    elif Rep1 > 0 and Rep2 < 0:
        type = "type_5_rep1_ins2"
    elif Rep1 == 0 and Rep2 > 0:
        type = "type_5_ref1_rep2"
    elif Rep1 == 0 and Rep2 == 0:
        type = "type_5_ref1_ref2"
    elif Rep1 == 0 and Rep2 > 0:
        type = "type_5_ref1_ins2"
    elif Rep1 < 0 and Rep2 > 0:
        type = "type_5_ins1_rep2"
    elif Rep1 < 0 and Rep2 == 0:
        type = "type_5_ins1_ref2"
    elif Rep1 < 0 and Rep2 < 0:
        type = "type_5_ins1_ins2"
    else:
        type = "other"
    if Rep3 > 0:
        type = type + "_rep3"
    elif Rep3 == 0:
        type = type + "_ref3"
    elif Rep3 < 0:
        type = type + "_ins3"
    else:
        type = "other"
    if Rep4 > 0:
        type = type + "_rep4"
    elif Rep4 == 0:
        type = type + "_ref4"
    elif Rep4 < 0:
        type = type + "_ins4"
    else:
        type = "other"

    new_line = "\t".join((query_id, length, type, str(Rep1), str(Rep2), str(Rep3), str(Rep4), strand1, query_start1, query_end1, sbjct_start1, sbjct_end1, strand2, query_start2, query_end2, sbjct_start2, sbjct_end2, strand3, query_start3, query_end3, sbjct_start3, sbjct_end3, strand4, query_start4, query_end4, sbjct_start4, sbjct_end4, strand5, query_start5, query_end5, sbjct_start5, sbjct_end5))
    print(new_line)

