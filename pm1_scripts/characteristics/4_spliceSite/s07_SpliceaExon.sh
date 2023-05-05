#!/bin/bash

fasta=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/genome/pm1.soft.fasta
#intronbed=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/annotation/pm1_scallop_multi-coding_intron.bed
#pre=pm1_intron
intronbed=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/annotation/pm1_scallop_final_intron.bed
pre=pm1_intron_full

cat $intronbed | awk -v OFS='\t' '{print $1,$2,$2+2,$4,$5,$6}' > ${pre}.rightmotif.bed
cat $intronbed | awk -v OFS='\t' '{print $1,$3-2,$3,$4,$5,$6}' > ${pre}.leftmotif.bed
        
wc -l $intronbed
wc -l ${pre}.leftmotif.bed
wc -l ${pre}.rightmotif.bed
	
#seperate Seq
bedtools getfasta -tab -fi $fasta -bed ${pre}.leftmotif.bed -fo ${pre}.leftmotif.txt
bedtools getfasta -tab -fi $fasta -bed ${pre}.rightmotif.bed -fo ${pre}.rightmotif.txt

#merge motif
cut -f 2 ${pre}.leftmotif.txt > ${pre}.leftmotif.tmp
cut -f 2 ${pre}.rightmotif.txt > ${pre}.rightmotif.tmp
paste -d '' ${pre}.leftmotif.tmp ${pre}.rightmotif.tmp > ${pre}.motif.seq
	
paste $intronbed ${pre}.leftmotif.txt ${pre}.rightmotif.txt ${pre}.motif.seq > ${pre}.motif.txt
sed -i '1i chr\tstart\tend\ttranscriptID\tScore\tstrand\tleftPosition\tleftSeq\trightPosition\trightSeq\tmotif' ${pre}.motif.txt


awk 'function rev(x) {r="";for(i=length(x);i;i--) r=r substr(x,i,1);return r} BEGIN{OFS="\t"} { if($6=="-") $11=rev($11); print $0}' ${pre}.motif.txt > ${pre}.motifRev.tmp
awk '{ if ($6 == "-") print }' ${pre}.motifRev.tmp | tr 'ACGT' 'TGCA' > ${pre}.motifRev.txt
awk '{ if ($6 == "+") print }' ${pre}.motifRev.tmp >> ${pre}.motifRev.txt
sed -i '1i chr\tstart\tend\ttranscriptID\tScore\tstrand\tleftPosition\tleftSeq\trightPosition\trightSeq\tmotifRev' ${pre}.motifRev.txt

rm ${pre}.leftmotif.bed ${pre}.leftmotif.txt ${pre}.rightmotif.bed ${pre}.rightmotif.txt ${pre}.motif.seq
rm ${pre}.leftmotif.tmp ${pre}.rightmotif.tmp
#rm ${pre}.motifRev.tmp

cat ${pre}.motifRev.txt | tr '[:lower:]' '[:upper:]' > ${pre}.motifRev.tsv
tail -n+2 ${pre}.motifRev.txt | awk '{print $11}' | tr '[:lower:]' '[:upper:]' | sort | uniq -c | tr -s " " | awk -v OFS='\t' '{print $2,$1,$1/119791}' | sort -k2 -n > ${pre}.motifRev.table
cp ${pre}.motifRev.table ../
awk '{sum+=$3} END {print sum}' ${pre}.motifRev.table
