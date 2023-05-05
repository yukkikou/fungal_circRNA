#!/bin/bash
#SBATCH -A ylab
#SBATCH -J rRNA
#SBATCH -p compute
#SBATCH -D /media/work/hxy/rRNA
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -o /home/hxy/log/rRNA-%A_%a.out
#SBATCH -e /home/hxy/log/rRNA-%A_%a.out
#SBATCH --mem 120G

workdir=/media/work/hxy/rRNA
fqdir=/media/share/node13/disk3/ylab/hxy/circFL_old/1_fq
rRNA=/media/share/node13/disk3/ylab/hxy/circFL_old/genome/rRNA/pm1rRNAnonETS
outdir=/media/share/node13/disk3/ylab/hxy/circFL_old/1_fq/1_rRNA

file=$fqdir/Y2.pass.fastq
thread=$SLURM_CPUS_PER_TASK

mkdir -p $workdir

bowtie2 -p $thread -x $rRNA \
	-r $file --no-unal \
	-S $workdir/$(basename $file ".pass.fastq").rRNA.sam &> $workdir/$(basename $file ".pass.fastq").rRNA.log

cp $workdir/$(basename $file ".pass.fastq").rRNA.log $outdir/
rm $workdir/$(basename $file ".pass.fastq").rRNA.sam  $workdir/$(basename $file ".pass.fastq").rRNA.log
