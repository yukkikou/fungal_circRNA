#!/uar/bin/python3

from collections import defaultdict
from collections import Counter
import re

def exonsplit(exons):
    exon_lst = exons.split(',')
    exon_ss = []
    exon_es = []
    for ex in exon_lst:
        exon_range = ex.split('|')[0]
        exon_ss.append(exon_range.split('-')[0])
        exon_es.append(exon_range.split('-')[1])
    return [exon_ss, exon_es]

def readlst(filename: str):
    with open(filename) as fin:
        for line in fin:
            cnts = line.rstrip().split('\t')
            circ_id = cnts[1]
            strand = cnts[3]
            exon_ss = exonsplit(cnts[4])[0]
            exon_es = exonsplit(cnts[4])[1]
            yield [circ_id, strand, exon_ss, exon_es]

def groupBycirc(lstIter):
    circs = defaultdict(list)
    for lstCnts in lstIter:
        try:
            circs[lstCnts[0]].append(lstCnts)
        except:
            raise ValueError(lstCnts)
    return circs

def exonSelect(circs, filename):
    for cId in circs:
        read_number = len(circs[cId])
        exon_number = []
        for i in range(read_number):
            exon_number.append(len(circs[cId][i][2]))
        block_number = int(Counter(exon_number).most_common(1)[0][0])
        block_reads = int(Counter(exon_number).most_common(1)[0][1])
        match_read = []
        for i in range(read_number):
            block_len = len(circs[cId][i][2])
            if block_len == block_number:
                match_read.append(circs[cId][i])
        exon_ss_pos = []
        exon_es_pos = []
        exon_size = []
        new_id = []
        for j in range(block_number):
            exon_ss_tmp = []
            exon_es_tmp = []
            exon_size_tmp = []
            for k in range(block_reads):
                exon_ss_tmp.append(match_read[k][2][j])
                exon_es_tmp.append(match_read[k][3][j])
            ss_keep_site = int(Counter(exon_ss_tmp).most_common(1)[0][0])-1
            es_keep_site = int(Counter(exon_es_tmp).most_common(1)[0][0])
            exon_size_tmp = str(es_keep_site-ss_keep_site)
            #print(exon_size_tmp)
            exon_size.append(exon_size_tmp)
            exon_ss_pos.append(str(ss_keep_site))
            exon_es_pos.append(str(es_keep_site))
            new_id = str(int(exon_ss_pos[0])+1)+'-'+str(exon_es_pos[-1])
        yield (match_read[0][0], block_number, ','.join(exon_size), ','.join(exon_ss_pos), ','.join(exon_es_pos), new_id)

if __name__ == '__main__':
    testfile = 'circlong_suppIN.lst'
    with open('circlong_suppIN.outv2', 'w') as fout:
        for cnts in exonSelect(groupBycirc(readlst(testfile))):
            print('\t'.join(str(x) for x in cnts), file = fout)
