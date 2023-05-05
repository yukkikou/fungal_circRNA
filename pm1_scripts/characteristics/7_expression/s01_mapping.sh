#!/usr/bin/bash

fqdir=/media/data7/hxy/PM1_illumima_202101/ANNO_XS01KF2020120029_PM-XS01KF2020120029-02_AHF3JTCCX2_2021-01-08/Rawdata/1_circIdentification
genome=/media/data5/hxy/genome/pm1.soft.fasta
hisatidx=/media/data5/hxy/genome/hisatindex/pm1
gtf=/media/data5/hxy/annotation/pm1_scallop_final.gtf

workdir=/media/data5/hxy/PM1_linear_illumina_2021
outdir=/media/data5/hxy/PM1_linear_illumina_2021
sifDir=/home/hxy/singularity

cldir=$workdir/1_clean
bamdir=$workdir/2_bam
fctsdir=$workdir/3_featureCounts

mkdir -p $workdir/{1_clean,2_bam,3_featureCounts}

thread=60
export SINGULARITY_BIND="/media:/media,/home/hxy:/home/hxy"

trimGaloreSif="singularity exec $sifDir/trimGalore.sif"
hisat2Sif="singularity exec $sifDir/hisat2.sif"
samtoolsSif="singularity exec $sifDir/samtools.sif"
subreadSif="singularity exec $sifDir/subread.sif"

#for sap in {Y2,Y3,EY3,EY4,C1,C2,EC1,EC3}
#do
#	raw1=$fqdir/$sap/1_fq/${sap}_R1.fq.gz
#	raw2=$fqdir/$sap/1_fq/${sap}_R2.fq.gz
#
#	cl1=$cldir/${sap}_R1_val_1.fq.gz
#	cl2=$cldir/${sap}_R2_val_2.fq.gz
#
#	bam=$bamdir/${sap}.hisat2.sorted.bam
#	tlog=$cldir/${sap}.trim_galore.log
#	ss=$bamdir/${sap}.hisat2.splice.site
#	hlog=$bamdir/${sap}.hisat2.summary
#
#	echo "*** galore start ***" && $trimGaloreSif trim_galore -j $thread -q 20 \
#    		--fastqc --fastqc_args "--nogroup" \
#		--gzip --length 20 \
#		-o $cldir \
#		--paired \
#		$raw1 $raw2 &> $tlog
#
#	echo "*** $sap: hisat2 start ***" && $hisat2Sif hisat2 -p $thread \
#	    --sensitive \
#	    --max-intronlen 2000 \
#	    --novel-splicesite-outfile $ss \
#	    --rf \
#	    -x $hisatidx \
#	    --summary-file $hlog \
#	    -1 $cl1 -2 $cl2 \
#        	| $samtoolsSif samtools sort -m 100G -o $bam -
#done

$subreadSif featureCounts -p -C -J \
    -T $thread -s 2 \
    -a $gtf -o $fctsdir/gene_count.tsv \
    -g 'gene_id' \
    $bamdir/*.bam

$subreadSif featureCounts -p -C -J \
    -T $thread -s 2 \
    -a $gtf -o $fctsdir/transcript_count.tsv \
    -g 'transcript_id' \
    $bamdir/*.bam

#mkdir -p $fctsdir/both
#$subreadSif featureCounts -p -C -B -J \
#    -T $thread -s 2 \
#    -a $gtf -o $fctsdir/both/gene_count_bothend.tsv \
#    -g 'gene_id' \
#    $bamdir/*.bam
#
#$subreadSif featureCounts -p -C -B -J \
#    -T $thread -s 2 \
#    -a $gtf -o $fctsdir/both/transcript_count_bothend.tsv \
#    -g 'transcript_id' \
#    $bamdir/*.bam
