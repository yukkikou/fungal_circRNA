#!/usr/bin/bash

circbed=IN_Quant.clean.bed
repbed=pm1.fasta.out.bed6
pre=flk

#win=500
# no reverse but not matter
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

genome=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/genome/pm1.soft.fasta
upflank=flk.up.bed
downflank=flk.down.bed

export SINGULARITY_BIND='/media:/media,/home/hxy:/home/hxy'

echo "*** Upfa ***" && singularity exec ~/singularity/bedtools/bedtools-2.30.0.sif \
    bedtools getfasta -s -nameOnly \
    -fi $genome \
    -bed $upflank \
    -fo | fold > $(basename $upflank ".bed").fa
echo "*** Downfa ***" && singularity exec ~/singularity/bedtools/bedtools-2.30.0.sif \
    bedtools getfasta -s -nameOnly \
    -fi $genome \
    -bed $downflank \
    -fo | fold > $(basename $downflank ".bed").fa


echo "*** MakeDB ***" && singularity exec ~/singularity/blast-2.13.0/blast-2.13.0.sif \
    makeblastdb -in $(basename $downflank ".bed").fa \
    -dbtype nucl -parse_seqids \
    -out downFlk

echo "*** Blastn ***" && singularity exec ~/singularity/blast-2.13.0/blast-2.13.0.sif \
    blastn -query $(basename $upflank ".bed").fa \
    -db downFlk \
    -out $(basename $upflank ".bed").rls \
    -evalue 1e-3 \
    -num_alignments 10 \
    -num_threads 20 \
    -outfmt "6 qseqid qstart qend sseqid sstart send evalue score pident nident qseq sseq"

echo "*** MakeDB ***" && singularity exec ~/singularity/blast-2.13.0/blast-2.13.0.sif \
    makeblastdb -in $(basename $upflank ".bed").fa \
    -dbtype nucl -parse_seqids \
    -out upFlk

echo "*** Blastn ***" && singularity exec ~/singularity/blast-2.13.0/blast-2.13.0.sif \
    blastn -query $(basename $downflank ".bed").fa \
    -db upFlk \
    -out $(basename $downflank ".bed").rls \
    -evalue 1e-3 \
    -num_alignments 10 \
    -num_threads 20 \
    -outfmt "6 qseqid qstart qend sseqid sstart send evalue score pident nident qseq sseq"

cat flk.down.rls | awk '$1==$4{print $0}' > flk.paired.rlt
