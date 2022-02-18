import sys

fin = open(sys.argv[1], "rt")
lines = fin.readlines()
for line in lines:
    fields = line.rstrip("\n").split("\t")
    query_id = fields.pop(0)
    length = fields.pop(0)
    type = fields.pop(0)
    count = len(fields)//5
    new_line = query_id + "\t" + length + "\t" + "type_1_Ref" + "\t" + "\t".join(fields)
    print(new_line)
