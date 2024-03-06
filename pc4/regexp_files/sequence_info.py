#!/usr/bin/python3
import re

def get_names(filename):
    with open(filename, "r") as file:
        file_con = file.read()

    return [seq.split(" ")[0] for seq in file_con.split(">")[1:]]

def get_most_common_OS(filename):
    with open(filename, "r") as file:
        file_con = file.read()

    dic = {}
    for os_name in re.finditer("OS=([^(OX=)]+)", file_con):
        name = file_con[os_name.span()[0]+3:os_name.span()[1]-1]
        if name in dic:
            dic[name] += 1
        else:
            dic[name] = 1

    return max(dic, key=dic.get)

def get_seq_dic(filename):
    with open(filename, "r") as file:
        file_con = file.read()
    
    seqs = file_con.split(">")[1:]
    seqs = [seq.split("\n")[0] for seq in seqs]

    dic = {}
    for seq in seqs:
        name = seq.split(" ")[0]
        os = re.search("OS=[^=]+", seq).span()
        os_val = seq[os[0]+3:os[1]-3]

        seq = seq[os[1]-2:]

        ox = re.search("OX=[^=]+", seq).span()
        ox_val = seq[ox[0]+3:ox[1]-3]

        seq = seq[ox[1]-2:]

        gn = re.search("GN=[^=]+", seq).span()
        gn_val = seq[gn[0]+3:gn[1]-3]

        seq = seq[gn[1]-2:]

        pe = re.search("PE=[^=]+", seq).span()
        pe_val = seq[pe[0]+3:pe[1]-3]

        sv_val = seq[pe[1]+1:]

        dic[name] = (os_val, ox_val, gn_val, pe_val, sv_val)
    
    return dic  

def get_pubmed(filename):
    with open(filename, "r") as file:
        file_con = file.read()
    
    return [file_con[pm.span()[0]+6:pm.span()[1]].strip() for pm in re.finditer("PUBMED(\s)+(\d)+", file_con)]

def clean_gb(filename):
    with open(filename, "r") as file:
        file_con = file.read()[6:]
    
    return re.sub("[\s\d]+", "", file_con).upper()


if __name__ == "__main__":
    # print(get_names("PS00727.fasta"))
    # print(get_most_common_OS("PS00727.fasta"))
    # dic = get_seq_dic("PS00727.fasta")
    # for key in dic:
    #     print("{:<21} {}".format(key, dic[key]))
    # print(get_pubmed("sequence.gb"))
    print(clean_gb("test.txt"))
