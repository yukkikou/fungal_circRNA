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
            chr_id = circ_id.split(':')[0]
            circ_ss = str(int(circ_id.split(':')[1].split('-')[0])-1)
            circ_es = str(circ_id.split(':')[1].split('-')[1])
            exon_ss = exonsplit(cnts[4])[0]
            exon_es = exonsplit(cnts[4])[1]
            yield [chr_id, circ_ss, circ_es, circ_id, strand, exon_ss, exon_es]

def groupBycirc(lstIter):
    circs = defaultdict(list)
    for lstCnts in lstIter:
        try:
            circs[lstCnts[3]].append(lstCnts)
        except:
            raise ValueError(lstCnts)
    return circs

def exonSelect(circs):
    for cId in circs:
        read_number = len(circs[cId])
        strand_tmp = []
        exon_number = []
        for i in range(read_number):
            strand_tmp.append(circs[cId][i][4].replace('NA', ''))
            exon_number.append(len(circs[cId][i][5]))
        strand_all = [x for x in strand_tmp if x != '']
        if len(strand_all) != 0:
            strand_final = Counter(strand_all).most_common(1)[0][0]
        else:
            strand_final = '.'
        block_number = int(Counter(exon_number).most_common(1)[0][0])
        block_reads = int(Counter(exon_number).most_common(1)[0][1])
        match_read = []
        for i in range(read_number):
            block_len = len(circs[cId][i][5])
            if block_len == block_number:
                match_read.append(circs[cId][i])
        for j in range(block_number):
            exon_ss_tmp = []
            exon_es_tmp = []
            for k in range(block_reads):
                exon_ss_tmp.append(match_read[k][5][j])
                exon_es_tmp.append(match_read[k][6][j])
            ss_keep_site = int(Counter(exon_ss_tmp).most_common(1)[0][0])-1
            es_keep_site = int(Counter(exon_es_tmp).most_common(1)[0][0])
            exon_size_tmp = str(es_keep_site-ss_keep_site)
            #print(exon_size_tmp)\
            yield (match_read[0][0], match_read[0][1], match_read[0][2], match_read[0][3], strand_final, block_number, exon_size_tmp, ss_keep_site, es_keep_site)

if __name__ == '__main__':
    testfile = '../0_cleanCircRNA/circlong_suppIN.lst'
    with open('circlong_suppIN.gtf.bed13', 'w') as fout:
        for cnts in exonSelect(groupBycirc(readlst(testfile))):
            print('\t'.join(str(x) for x in cnts), file = fout)
