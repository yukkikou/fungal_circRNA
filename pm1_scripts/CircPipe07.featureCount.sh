#!/usr/bin/bash
#SBATCH -p single-ib
#SBATCH --node 1
#SBATCH --cpu-per-tasks=20
#SBATCH -J fcts
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%a_%a.out
#SBATCH --mem 120G
#SBATCH --array 1

gtf=/media/share/node10/disk3/hxy/reference/Talaromyces_marneffei/annotation/pm1_scallop_final.gtf
workdir=/media/work/hxy/PM1_linear_illumina_2021/6_linearQuantify
outdir=/media/share/node10/disk3/ylab/hxy/PM1_linear_illumina_2021
bamdir=$outdir/2_bam
bams=`find $bamdir -name "*.sorted.bam"`

sifDir=/media/share/node09/disk1/ylab/dmh/BY2022
thread=$SLURM_CPUS_PER_TASK

fctsSif="singularity exec $sifDir/subread.sif"
export SINGULARITY_BIND="/media:/media,/home/hxy:/home/hxy"

# gene
$fcts featureCounts -p -B -C \
    -T $thread -s 2 \
    -a $gtf -o $workdir/gene_count.tsv \
    -g 'gene_id' \
    $bams

# transcript
$fcts featureCounts -p -B -C \
    -T $thread -s 2 \
    -a $gtf -o $workdir/transcript_count.tsv \
    -g 'transcript_id' \
    $bams
 

