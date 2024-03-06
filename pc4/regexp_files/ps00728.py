#!/usr/bin/python3
import re
import sys

class Seq:
    def __init__(self, seq):
        seq = seq.split("\n")
        self.name = seq[0]
        self.seq = "".join(seq[1:])

def load_sequences(filename):
    ret_seqs = []
    with open(filename, "r") as file:
        file_seqs = file.read()
        seqs = file_seqs.split(">")
        for seq in seqs:
            if (len(seq) > 0):
                ret_seqs.append(Seq(seq))
            
    return ret_seqs

def construct_motif():
    return "N.G.R[LIVM]D[LIVMFYH].[LV].S"


def find_motif(seq, motif):
    import prosite_re as pr
    mo = re.search(motif, seq.seq)

    if (mo != None):
        return mo.span()
    else:
        return None
    
def detect_motifs(seqs):
    for seq in seqs:
        motif_span = find_motif(seq, construct_motif())
        
        if motif_span:
            print("\033[92m{:<21} MATCH\033[0m".format(seq.name.split(" ")[0]))
        else:
            print("\033[91m{:<21} NOT_MATCH\033[0m".format(seq.name.split(" ")[0]))


if __name__ == '__main__':
    seqs = load_sequences(sys.argv[1])
    detect_motifs(seqs)
    


