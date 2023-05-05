#!/usr/bin/bash

circbed=IN_Quant.clean.bed
repbed=pm1.fasta.repeat.out.bed6
pre=flk

#win=500

cat $circbed | awk 'BEGIN{OFS="\t";} {
	if ($6=="+"){
		print $1,$2-501,$2,$4,$5,$6;
	} else {
		print $1,$3,$3+500,$4,$5,$6;
	}
}' > ${pre}.up.bed

cat $circbed | awk 'BEGIN{OFS="\t";} {
	if ($6=="+"){
		print $1,$3,$3+500,$4,$5,$6
	} else {
		print $1,$2-501,$2,$4,$5,$6
	}
}' > ${pre}.down.bed


bedtools intersect -a ${pre}.up.bed -b $repbed -wo > ${pre}.up.rep.inter
bedtools intersect -a ${pre}.down.bed -b $repbed -wo > ${pre}.down.rep.inter


