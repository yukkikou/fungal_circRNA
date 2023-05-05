#!/usr/bin/env bash
#SBATCH -p single-ib,single-em
#SBATCH --nodes 1
#SBATCH --cpus-per-task=35
#SBATCH -J hisat2
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.err
#SBATCH --mem 120G
#SBATCH --array 0-14:2

samlst=/media/share/node10/disk3/ylab/hxy/PM1_circSeq/PM1_illumima_202101/circ_sample.list
fqdir=/media/share/node11/disk3/ylab/hxy/PM1_circSeq/PM1_illumima_202101
genome=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/genome/pm1.soft.fasta
hisatidx=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/index/PM1
gtf=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/annotation/pm1_scallop_final.gtf

workdir=/media/work/hxy/PM1_circ_nanopore_202101/6_linearQuantify
outdir=/media/share/node10/disk3/ylab/hxy/PM1_linear_illumina_2021
sifDir=/media/share/node09/disk1/ylab/dmh/BY2022

cldir=$workdir/fq
bamdir=$workdir/bam
tmpdir=$workdir/tmp
mkdir -p ${workdir}/{fq,bam,tmp}

thread=$SLURM_CPUS_PER_TASK
nmem=$SBATCH_MEM_PER_NODE

i=$SLURM_ARRAY_TASK_ID
circ_array=(`cat $samlst`)
sample_array=(${circ_array[@]:$i:2})
raw_prefix=${sample_array[0]}
clean_prefix=${sample_array[0]}

raw1=$fqdir/${raw_prefix}_R1.fq.gz
raw2=$fqdir/${raw_prefix}_R2.fq.gz

cl1=$cldir/${clean_prefix}_R1_val_1.fq.gz
cl2=$cldir/${clean_prefix}_R2_val_2.fq.gz

bam=$bamdir/${clean_prefix}.hisat2.sorted.bam
tlog=$cldir/${clean_prefix}.trim_galore.log
ss=$bamdir/${clean_prefix}.hisat2.splice.site
hlog=$bamdir/${clean_prefix}.hisat2.summary

export SINGULARITY_BIND="/media:/media,/home/hxy:/home/hxy"

trimGaloreSif="singularity exec $sifDir/trimGalore.sif"
hisat2Sif="singularity exec $sifDir/hisat2.sif"
samtoolsSif="singularity exec $sifDir/samtools.sif"
subreadSif="singularity exec $sifDir/subread.sif"

echo "*** galore start ***" && $trimGaloreSif trim_galore -j $thread -q 20 \
    --fastqc --fastqc_args "--nogroup" \
    --gzip --length 20 \
    -o $cldir \
    --paired \
    $raw1 $raw2 &> $tlog

echo "*** hisat2 start ***" && $hisat2Sif hisat2 -p $thread \
    --sensitive \
    --max-intronlen 2000 \
    --novel-splicesite-outfile $ss \
    --rf \
    -x $hisatidx \
    --summary-file $hlog \
    -1 $cl1 -2 $cl2 \
        | $samtoolsSif samtools sort -m 40G -o $bam -

# data transfor
mkdir -p $outdir/1_clean/${clean_prefix}
mkdir -p $outdir/2_bam/${clean_prefix}

cp $cldir/${clean_prefix}* $outdir/1_clean/${clean_prefix}/
cp $bam $ss $hlog $outdir/2_bam/${clean_prefix}/

# work space release
rm $cldir/${clean_prefix}* $bam $ss $hlog 
