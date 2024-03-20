# -*- coding: utf-8 -*-

class MyBlast:
    '''
    Classe para matrizes de pontos
    '''

    def __init__(self, filename = None, w = 3):
        '''
        Construtor
        '''
        if filename is not None:
            self.readDatabase(filename)
        else:
            self.db = []
        self.w = w
        self.map = None

    def readDatabase(self, filename):
        """From file with sequences line by line read the sequences to a list"""
        lines = []
        with open(filename, "r") as f:
            lines = f.readlines()

        self.db = lines

    def setDatabase(self, db):
        self.db = db

    def addSequenceDB(self, seq):
        """Add an extra sequence to DB"""
        self.db.append(seq)

    def buildMap (self, query):
        res = {}
        for i in range(len(query) - self.w + 1):
            seq = query[i:i + self.w]
            if seq in res:
                res[seq].append(i)
            else:
                res[seq] = [i]
                
        self.map = res

    def getHits (self, seq, query): 
        res = []
        for i in range(len(seq) - self.w + 1):
            sub_seq = seq[i:i + self.w]
            if sub_seq in self.map:
                for index in self.map[sub_seq]:
                    res.append((index, i))

        return res
    
    def getHitsMiss (self, seq, query, missmatch): 
        res = []
        for i in range(len(seq) - self.w + 1):
            sub_seq = seq[i:i + self.w]
            for map_seq in self.map.keys():
                if hamming_distace(map_seq, sub_seq) < missmatch:
                    for index in self.map[sub_seq]:
                        res.append((index, i))

        return res

    def extendsHit (self, seq, hit, query):
        stq, sts = hit[0], hit[1]
        matfw = 0
        k=0
        bestk = 0
        while 2*matfw >= k and stq+self.w+k < len(query) and sts+self.w+k < len(seq):
            if query[stq+self.w+k] == seq[sts+self.w+k]:
                matfw+=1
                bestk = k+1
            k += 1
        size = self.w + bestk

        k = 0
        matbw = 0
        bestk = 0
        while 2*matbw >= k and stq > k and sts > k:
            if query[stq-k-1] == seq[sts-k-1]:
                matbw+=1
                bestk = k+1
            k+=1
        size += bestk

        return (stq-bestk, sts-bestk, size, self.w+matfw+matbw)

    def hitBestScore(self, seq, query):
        hits = self.getHits(seq, query)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = self.extendsHit(seq, h, query)
            score = ext[3]
            if score > bestScore or (score== bestScore and ext[2] < best[2]):
                bestScore = score
                best = ext
        return best

    def bestAlignment (self, query):
        self.buildMap(query)
        bestScore = -1.0
        res = (0,0,0,0,0)
        for k in range(0,len(self.db)):
            bestSeq = self.hitBestScore(self.db[k], query)
            if bestSeq != ():
                score = bestSeq[3]
                if score > bestScore or (score== bestScore and bestSeq[2] < res[2]):
                    bestScore = score
                    res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
        if bestScore < 0: return ()
        else: return res

def test1():
    mb = MyBlast("./sequences/seqBlast.txt", 11)
    query = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
    r = mb.bestAlignment(query)
    print(r)

def test2():
    mb = MyBlast("./sequences/seqBlast.txt", 11)
    query2 = "cgacgacgacgacgaatgatg"
    r = mb.bestAlignment(query2)
    print(r)
    
def extractDB(filename):
    res = []
    with open(filename) as f:
        dbs = f.read().split(">")[1:]
        for db in dbs:
            res_db = db.split("\n")[1:]
            res.append("".join(res_db))
    return res


def test3():
    import sequence_alignments as sa

    mb = MyBlast(None, 11)
    mb.setDatabase(extractDB("./glyco_sequences/db.fas"))
    query = extractDB("./glyco_sequences/query.fas")[0]
    r = mb.bestAlignment(query)
    print(r)
    # query3 = MyBlast("./glyco_sequences/query.fas", 11)
    sm =sa.create_submat(2, -1, "ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    seq1 = query
    seq2 = mb.db[r[4]]
    S, T = sa.needleman_Wunsch(seq1, seq2, sm, -8)
    print(sa.recover_align(T, seq1, seq2))
    # r2 = sa.smith_Waterman(mb, query3)
    # print(r1,r2)

test3()