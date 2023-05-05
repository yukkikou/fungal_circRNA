#!/usr/bin/bash

genome=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/genome/pm1.soft.fasta
#circbed=IN_Quant.clean.bed
exonbed=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/annotation/pm1_scallop_final_exon.bed
#pre=flk
pre=exon
#win=20

export SINGULARITY_BIND='/media:/media,/home/hxy:/home/hxy'

cat $exonbed | awk -v OFS='\t' '$2>22{print $0}' > pm1_scallop_final_exon.clean.bed
circbed=pm1_scallop_final_exon.clean.bed

# for minus strand, reverse complementary are needed
cat $circbed | awk 'BEGIN{OFS="\t";} {
    if ($6=="+"){
        print $1,$2-22,$2+20,$4,45,$6
    }
}' > ${pre}.positive.up.bed

cat $circbed | awk 'BEGIN{OFS="\t";} {
    if ($6=="+"){
        print $1,$3-20,$3+22,$4,$5,$6
    }
}' > ${pre}.positive.down.bed

cat $circbed | awk 'BEGIN{OFS="\t";} {
    if ($6=="-"){
        print $1,$3-20,$3+22,$4,$5,$6
    }
}' > ${pre}.minus.up.bed

cat $circbed | awk 'BEGIN{OFS="\t";} {
    if ($6=="-"){
        print $1,$2-22,$2+20,$4,$5,$6
    }
}' > ${pre}.minus.down.bed


echo "*** Positive_Up_fasta ***" && singularity exec ~/singularity/bedtools/bedtools-2.30.0.sif \
    bedtools getfasta -s -nameOnly \
    -fi $genome \
    -bed ${pre}.positive.up.bed \
    -fo | fold > ${pre}.positive.up.fa

echo "*** Positive_down_fasta ***" && singularity exec ~/singularity/bedtools/bedtools-2.30.0.sif \
    bedtools getfasta -s -nameOnly \
    -fi $genome \
    -bed ${pre}.positive.down.bed \
    -fo | fold > ${pre}.positive.down.fa

echo "*** Minus_up_fasta ***" && singularity exec ~/singularity/bedtools/bedtools-2.30.0.sif \
    bedtools getfasta -nameOnly \
    -fi $genome \
    -bed ${pre}.minus.up.bed \
    -fo | fold > ${pre}.minus.up.fa.raw

singularity exec ~/singularity/seqkit/seqkit-2.2.0.sif seqkit seq \
    -p -u -r -t dna -j 10 \
    ${pre}.minus.up.fa.raw > ${pre}.minus.up.fa

echo "*** Minus_down_fasta ***" && singularity exec ~/singularity/bedtools/bedtools-2.30.0.sif \
    bedtools getfasta -nameOnly \
    -fi $genome \
    -bed ${pre}.minus.down.bed \
    -fo | fold > ${pre}.minus.down.fa.raw

singularity exec ~/singularity/seqkit/seqkit-2.2.0.sif seqkit seq \
    -p -u -r -t dna -j 10 \
    ${pre}.minus.down.fa.raw > ${pre}.minus.down.fa

cat ${pre}.positive.up.fa ${pre}.minus.up.fa > ${pre}.up.fa
cat ${pre}.positive.down.fa ${pre}.minus.down.fa > ${pre}.down.fa
rm ${pre}.minus.up.fa.raw ${pre}.minus.down.fa.raw
