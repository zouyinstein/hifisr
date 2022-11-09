import sys

sample=sys.argv[1]
with open("all_" + sys.argv[2] + ".ids", "rt") as fin:
    ids = fin.readlines()

for id in ids:
    ID = "__".join(id.rstrip("\n").split("/"))
    with open(sys.argv[2] + "_commands/runme_" + ID + ".sh", "at") as fout:
        print("blastn -query " + sample + ".sorted.REF_" + sys.argv[2] + ".fasta.split/" + sample + ".sorted.REF_" + sys.argv[2] + ".id_" + ID + ".fasta -db " + sample + "_" + sys.argv[2] + ".fa -outfmt 5 -out tmp/blastn_" + sample + "_" + sys.argv[2] + ID + "_" + ".xml", file=fout)
        print("python ../../scripts/analysis_blastn.py tmp/blastn_" + sample + "_" + sys.argv[2] + "_" + ID + "_" + ".xml " + sample + "_" + sys.argv[2] + ".fa > tmp/blastn_" + sample + "_" + sys.argv[2] + ID + ".txt", file=fout)
