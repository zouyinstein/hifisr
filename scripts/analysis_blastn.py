from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys

class Blastn_info():
    def __init__(self, query_id, length, alignments):
        self.query_id = query_id
        self.length = length
        self.alignments = alignments

    def merge_border_hits(self, ref_length):
        length = len(self.alignments)
        for i in range(0, len(self.alignments)):
            for j in range(0, len(self.alignments)):
                if (self.alignments[i][2] + 1) == self.alignments[j][1]:
                    strand1 = self.alignments[i][0]
                    strand2 = self.alignments[j][0]
                    if strand1 == "Plus" and strand2 == "Plus" and self.alignments[i][4] == ref_length and self.alignments[j][3] == 1:
                        merged_query_start = self.alignments[i][1]
                        merged_query_end = self.alignments[j][2]
                        merged_sbjct_start = self.alignments[i][3]
                        merged_sbjct_end = self.alignments[j][4]
                        if i < j:
                            self.alignments.pop(j)
                            self.alignments.pop(i)
                        else:
                            self.alignments.pop(i)
                            self.alignments.pop(j)
                        self.alignments.append(("Plus", merged_query_start, merged_query_end, merged_sbjct_start, merged_sbjct_end))
                        return
                    elif strand1 == "Minus" and strand2 == "Minus" and self.alignments[i][4] == 1 and self.alignments[j][3] == ref_length:
                        merged_query_start = self.alignments[i][1]
                        merged_query_end = self.alignments[j][2]
                        merged_sbjct_start = self.alignments[i][3]
                        merged_sbjct_end = self.alignments[j][4]
                        if i < j:
                            self.alignments.pop(j)
                            self.alignments.pop(i)
                        else:
                            self.alignments.pop(i)
                            self.alignments.pop(i)
                        self.alignments.append(("Minus", merged_query_start, merged_query_end, merged_sbjct_start, merged_sbjct_end))
                        return

    def get_piece(self):
        # pop one line per call
        strand = self.alignments[0][0]
        query_start = self.alignments[0][1]
        query_end = self.alignments[0][2]
        sbjct_start = self.alignments[0][3]
        sbjct_end = self.alignments[0][4]
        for i in range(1, len(self.alignments)):
            if self.alignments[i][1] <= self.alignments[0][1] and self.alignments[0][2] <= self.alignments[i][2]:
                self.alignments.pop(0)
                return "yes"
            elif self.alignments[0][1] <= self.alignments[i][1] and self.alignments[i][2] <= self.alignments[0][2]:
                self.alignments.pop(i)
                return "yes"
        self.alignments.pop(0)
        return (strand, query_start, query_end, sbjct_start, sbjct_end)
    
    def sort_pieces(self):
        return Blastn_info(self.query_id, self.length, self.alignments.sort(key=lambda x:x[1]))

    def write_pieces(self):
        list_alignment = list()
        for alignment in self.alignments:
            str_alignment = "\t".join((alignment[0], str(alignment[1]), str(alignment[2]), str(alignment[3]), str(alignment[4])))
            list_alignment.append(str_alignment)
        line = "\t".join((self.query_id, str(self.length), "\t".join(list_alignment)))
        print(line)
    # calculate repsize, gapsize : in another script
    # check the type: in another script


def get_blastn_info(blast_record): 
    query_id = blast_record.query
    length = blast_record.query_length  # type int
    for alignment in blast_record.alignments:
        # line_num = len(alignment.hsps)
        all_strand = list()
        all_query_start = list() # type int
        all_query_end = list()
        all_sbjct_start = list()
        all_sbjct_end = list()
        for hsp in alignment.hsps:
            all_strand.append(hsp.strand[1])  # 'Plus', 'Minus'
            all_query_start.append(hsp.query_start)
            all_query_end.append(hsp.query_end)
            all_sbjct_start.append(hsp.sbjct_start)
            all_sbjct_end.append(hsp.sbjct_end)
    alignments = list(zip(all_strand, all_query_start, all_query_end, all_sbjct_start, all_sbjct_end))
    return Blastn_info(query_id, length, alignments)


if __name__ == "__main__":
    result_handle = open(sys.argv[1])
    ref_length = len(SeqIO.read(sys.argv[2], "fasta").seq)
    if sys.argv[2].endswith("_plastid.fa"):
        # Currently, do not merge cross-border reads for platid genome because of the very large repeat IRb.
        ref_length = 1000000
    blast_records = list(NCBIXML.parse(result_handle))
    for blast_record in blast_records:
        # get info from one  record
        blastn_info = get_blastn_info(blast_record)
        # merge continuous lines across borders
        count = len(blastn_info.alignments)
        blastn_info.merge_border_hits(ref_length)
        while len(blastn_info.alignments) < count:
            count = len(blastn_info.alignments)
            blastn_info.merge_border_hits(ref_length)
        # get pieces from merged blastn_info
        blastn_pieces = Blastn_info(blastn_info.query_id, blastn_info.length, list())
        piece = "yes"
        while len(blastn_info.alignments) > 0:
            piece = blastn_info.get_piece()
            while piece == "yes":
                piece = blastn_info.get_piece()
            blastn_pieces.alignments.append(piece)
        # sort blastn_pieces
        blastn_pieces.sort_pieces()
        # print to screen
        blastn_pieces.write_pieces()