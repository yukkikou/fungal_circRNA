#!/bin/bash
#SBATCH -A ylab
#SBATCH -J col_ZFQ
#SBATCH -p compute
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.out
#SBATCH --mem 100G

# mamba env: CIRI-long
genome=/media/share/node10/disk3/ylab/hxy/reference/Mus_musculus/Genome/GRCm38.p4.genome.fa
gtf=/media/share/node10/disk3/ylab/hxy/reference/Mus_musculus/Annotation/gencode.vM10.annotation.gtf

workdir=/media/work/hxy/CIRIlong_zfq/6_cirilong/2_collapse
cldir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/5_nano
outdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/6_cirilong

circlst=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/script/circFA.lst
logdir=$workdir/log
prefix=Mus

mkdir -p $logdir

thread=$SLURM_CPUS_PER_TASK

# conda env: CIRI-long
CIRI-long collapse -i $circlst \
	-o $workdir \
	-r $genome \
	-p $prefix \
	-a $gtf \
	-t $thread &> $logdir/${prefix}.collapes.log

# data transfor
mkdir -p $outdir/log

cp -r $workdir $outdir
cp $logdir/${prefix}.collapes.log $outdir/log

# work space release
rm -r $workdir
rm $logdir/${prefix}.collapes.log
