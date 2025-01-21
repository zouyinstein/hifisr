import sys

with open(sys.argv[1], "rt") as fin:
    lines = [ line.rstrip("\n") for line in fin.readlines() ]

for line in lines:
    fields = line.split("\t")
    if len(fields) == 8:
        if (int(fields[5]) - int(fields[4]))/int(fields[1]) > 0.99:
            print(fields[0])
    else:
        with open(sys.argv[1].split(".")[0] + ".error", "at") as fout:
            print(fields[0], file=fout)



