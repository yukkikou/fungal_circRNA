#!/usr/bin/python3

from collections import defaultdict

def readReads(filename):
    with open(filename) as fin:
        for line in fin:
            cnts = line.rstrip().split('\t')
            cntSpli = cnts[4].split(',')
            #print(cntSpli)
            yield [cnts[1],cnts[3],cntSpli]

def splitLoci(readItr):
    loci = defaultdict(list)
    for circ in readItr:
        loci[circ[0]].append(circ)
        #print(loci[circ[0]][0])
    return loci

def writeBed(filename, loci):
    with open(filename, "w") as fout:
        for cLoci in loci:
            #print(cLoci)
            readNum = len(loci[cLoci])
            print('RN: ' + str(readNum))
            for i in range(readNum):
                featNum = len(loci[cLoci][i][2])
                for j in range(featNum):
                    #print(loci[cLoci][i][0], loci[cLoci][i][1], loci[cLoci][i][2][j], sep = '\t')
                    chrn = loci[cLoci][i][0].split(":")[0]
                    site = loci[cLoci][i][2][j].split("|")[0].split('-')
                    ssite = int(site[0])-1
                    esite = int(site[1])
                    cId = loci[cLoci][i][0]
                    score = chrn+'_'+loci[cLoci][i][2][j].split("|")[0]+'_'+str(featNum)
                    strand = loci[cLoci][i][1]
                    print(chrn, ssite, esite, cId, score, strand, sep = '\t', file = fout)

if __name__ == "__main__":
    testfile='circlong_suppIN.lst'
    writeBed('circlong_suppIN.bed',splitLoci(readReads(testfile)))
