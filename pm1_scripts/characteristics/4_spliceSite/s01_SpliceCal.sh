#!/bin/bash
# For illumina motif calculate

fasta=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/genome/pm1.soft.fasta
gtfdir=/media/share/node10/disk3/ylab/hxy/PM1_circ_nanopore_202101/5_BSJsite/1_quantify/3_All_l2

for pre in {IN_M1,IN_M2,IN_M3,IN_M4,IN_Y1,IN_Y2,IN_Y3,IN_Y4}
do
        tail -n +6 $gtfdir/$pre/${pre}_Quant.gtf | awk -v OFS='\t' '{print $1,$4,$5,$7,$10}' | sed 's/"//g;s/;//g'> ${pre}_Quant.bed
        tail -n +6 $gtfdir/$pre/${pre}_Quant.gtf | awk -v OFS='\t' '{print $1,$4-3,$4-1,$7,$10}'> ${pre}.leftmotif.bed
	tail -n +6 $gtfdir/$pre/${pre}_Quant.gtf | awk -v OFS='\t' '{print $1,$5,$5+2,$7,$10}'> ${pre}.rightmotif.bed
        echo '-------'${pre}'--------'
        wc -l ${pre}_Quant.bed
        wc -l ${pre}.leftmotif.bed
	wc -l ${pre}.rightmotif.bed
        #truely bed: 0 based
	
	#seperate Seq
	bedtools getfasta -tab -fi $fasta -bed ${pre}.leftmotif.bed -fo ${pre}.leftmotif.txt.raw
	bedtools getfasta -tab -fi $fasta -bed ${pre}.rightmotif.bed -fo ${pre}.rightmotif.txt.raw

    # to upper
    cat ${pre}.leftmotif.txt.raw | tr [:lower:] [:upper:] > ${pre}.leftmotif.txt
    cat ${pre}.rightmotif.txt.raw | tr [:lower:] [:upper:] > ${pre}.rightmotif.txt
    
    #merge motif
	cut -f 2 ${pre}.leftmotif.txt > ${pre}.leftmotif.tmp
	cut -f 2 ${pre}.rightmotif.txt > ${pre}.rightmotif.tmp
	paste -d '' ${pre}.leftmotif.tmp ${pre}.rightmotif.tmp > ${pre}.motif.seq
	
	paste ${pre}_Quant.bed ${pre}.leftmotif.txt ${pre}.rightmotif.txt ${pre}.motif.seq > ${pre}.motif.txt
	sed -i '1i chr\tstart\tend\tstrand\tcircID\tleftPosition\tleftSeq\trightPosition\trightSeq\tmotif' ${pre}.motif.txt

	awk 'function rev(x) {r="";for(i=length(x);i;i--) r=r substr(x,i,1);return r} BEGIN{OFS="\t"} { if($4=="-") $10=rev($10); print $0}' ${pre}.motif.txt > ${pre}.motifRev.tmp
	awk '{ if ($4 == "-") print }' ${pre}.motifRev.tmp | tr 'ACGT' 'TGCA' > ${pre}.motifRev.txt
	awk '{ if ($4 == "+") print }' ${pre}.motifRev.tmp >> ${pre}.motifRev.txt
	sed -i '1i chr\tstart\tend\tstrand\tcircID\tleftPosition\tleftSeq\trightPosition\trightSeq\tmotifRev' ${pre}.motifRev.txt

	rm ${pre}.leftmotif.bed ${pre}.leftmotif.txt ${pre}.rightmotif.bed ${pre}.rightmotif.txt ${pre}.motif.seq
	rm ${pre}.leftmotif.tmp ${pre}.rightmotif.tmp
	rm ${pre}.motifRev.tmp
done

cat *.motifRev.txt | sort | uniq > all.motifRev.txt
cat all.motifRev.txt | tr [:lower:] [:upper:] > all.motifRev.tsv
cat all.motifRev.txt | cut -f10 | sed '1d' | tr [:lower:] [:upper:] | sort | uniq -c > all.motifRev.table
