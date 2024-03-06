#!/usr/bin/python3

def read_fasta(filename):
    with open(filename, "r") as file:
        file_seqs = file.readlines()[1:]

    return "".join(file_seqs).replace("\n","")

def find_zync_finger(seq):
    import re
    spans = []
    regexp = "C.H.[LIVMFY]C.{2}C[LIVMYA]"
    for mo in re.finditer(regexp, seq):
        spans.append(mo.span())

    return spans

    
def find_prosite(seq, profile):
    from re import search
    regexp = profile.replace("-","")
    regexp = regexp.replace("x",".")
    regexp = regexp.replace("(","{")
    regexp = regexp.replace(")","}")
    mo = search(regexp, seq)
    
    if (mo != None):
        return mo.span()[0]
    else:
        return -1
    
def highlight_motif(seq, span):
    seq = seq.lower()
    for sp in span:
        seq = seq[0:sp[0]] +"\033[92m"+ seq[sp[0]:sp[1]].upper() + "\033[0m" + seq[sp[1]:]
    return seq
    
def test():
    seq = "HKMMLASCKHLLCLKCIVKLG"
    seq = read_fasta("Q8RXD4.fasta.txt")
    zync_span = find_zync_finger(seq)
    print(zync_span)
    print(find_prosite(seq,"C-x-H-x-[LIVMFY]-C-x(2)-C-[LIVMYA]"))
    print(highlight_motif(seq, zync_span))

if __name__ == '__main__':    
    test()