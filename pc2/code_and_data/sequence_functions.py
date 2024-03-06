#!/usr/bin/python3
# -*- coding: utf-8 -*-

def read_seq_from_file(filename):
    """ Reads a sequence from a multi-line text file. """
    fh = open(filename, "r")
    lines = fh.readlines()
    seq = ""
    for l in lines:
        seq += l.replace("\n","")
    fh.close()
    return seq

def write_seq_to_file(seq, filename):
    """ Writes a sequence to file. """
    file = open(filename, "w")
    file.write(seq)
    file.close()

def read_genetic_code_from_file(filename):
    """ Reads the genetic code to a dictionary from a multi-line text file. """
    dic = {}
    seq = read_seq_from_file(filename)[1:]
    lines = seq.split("\"\"")
    for line in lines:
        words = line.split(" ")
        codon = words[0][0:-1]
        aa = words[1][1:]
        dic[codon] = aa

    return dic

def validate_dna (dna_seq):
    """ Checks if DNA sequence is valid. Returns True is sequence is valid, or False otherwise. """
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")
    if valid == len(seqm): return True
    else: return False

def frequency (seq):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """
    dic = {}
    for s in seq.upper():
        if s in dic: dic[s] += 1
        else: dic[s] = 1
    return dic

def gc_content (dna_seq):
    """ Returns the percentage of G and C nucleotides in a DNA sequence. """
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc": gc_count += 1
    return gc_count / len(dna_seq)

def gc_content_subseq (dna_seq, k=100):
    """ Returns GC content of non-overlapping sub-sequences of size k. """
    # complete
    # ...


def transcription (dna_seq):
    """ Function that computes the RNA corresponding to the transcription of the DNA sequence provided. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    return dna_seq.upper().replace("T","U")


def reverse_complement (dna_seq):
    """ Computes the reverse complement of the inputted DNA sequence. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    # complete function
    # A > T
    # T > A
    # C > G
    # G > C
    dic = {"A":"T", "T":"A", "C":"G", "G":"C"}
    comp_rev = ""
    for s in dna_seq:
        if s in dic:
            comp_rev = dic[s] + comp_rev
        else:
            comp_rev = "X" + comp_rev

    return comp_rev

def translate_codon (cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc: return tc[cod]
    else: return None


def translate_seq (dna_seq, ini_pos = 0):
    """ Translates a DNA sequence into an aminoacid sequence. """
    # ini_pos = 0 > frame 1
    # ini_pos = 1 > frame 2
    # ....
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    seq_aa = ""
    for i in range(ini_pos, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        aa = translate_codon(cod)
        if aa: seq_aa += aa
        else: seq_aa += "X"

    return seq_aa

def codon_usage(dna_seq, aa):
    """Provides the frequency of each codon encoding a given aminoacid, in a DNA sequence ."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    dic = {}
    total = 0
    for i in range(0, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic:
                dic[cod] += 1
            else: dic[cod] = 1
            total += 1
    if total >0:
        for k in dic:
            dic[k] /= total
    return dic


def reading_frames (dna_seq):
    """Computes the six reading frames of a DNA sequence (including the reverse complement."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    res.append(translate_seq(dna_seq,0))
    res.append(translate_seq(dna_seq,1))
    res.append(translate_seq(dna_seq,2))
    rc = reverse_complement(dna_seq)
    res.append(translate_seq(rc,0))
    res.append(translate_seq(rc,1))
    res.append(translate_seq(rc,2))
    return res


def all_proteins_rf (aa_seq):
    """Computes all posible proteins in an aminoacid sequence."""
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins

def all_orfs (dna_seq):
    """Computes all possible proteins for all open reading frames."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    for rf in reading_frames(dna_seq):
        prots = all_proteins_rf(rf)
        for p in prots:
            res.append(p)
    return res

def all_orfs_ord (dna_seq, minsize = 0):
    """Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames (dna_seq)
    res = []
    orfs = all_orfs(dna_seq)
    for p in orfs:
        if len(p) >= minsize:
            insert_prot_ord(p, res)
    return res

def insert_prot_ord (prot, list_prots):
    ''' inserts prot in list_prots in a sorted way '''
    i = 0
    while i < len(list_prots) and len(prot) < len(list_prots[i]):
        i += 1
    list_prots.insert(i, prot)


def test_frequency():
    seq_aa = input("Protein sequence:")
    freq_aa = frequency(seq_aa)
    list_f = sorted(freq_aa.items(), key=lambda x: x[1], reverse = True)
    for (k,v) in list_f:
        print("Aminoacid:", k, ":", v)

def test_sequence_from_file(filename):
    seq = read_seq_from_file(filename)
    gc = gc_content(seq)
    ref = reading_frames(seq)
    return (gc, ref)

def read_fasta_2dictionary(filename):
    file = open(filename, "r")
    sequences = file.read().split(">")[1:]
    file.close()

    dic = {}

    for seq in sequences:
        lines = seq.split("\n")
        name = lines[0].split(" ")[0]

        sequence = ""
        for l in lines[1:]:
            sequence += l.replace("\n","")

        dic[name] = sequence

    return dic

if __name__ == "__main__":
    # test here all implemented functions
    # used your own defined sequences or read from example files

    ## add here test function 
    # print(test_sequence_from_file("example_Hinfluenzae.txt"))
    # print(all_orfs(read_seq_from_file("example_Hinfluenzae.txt")))
    # print(all_orfs_ord(read_seq_from_file("example_Hinfluenzae.txt")))
    # print(len(all_orfs_ord(read_seq_from_file("example_Hinfluenzae.txt"))[0]))
    # print(read_fasta_2dictionary("PS00727.fasta"))
    l = "abcd"

    print(l[::-1])
