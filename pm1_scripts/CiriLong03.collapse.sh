#!/bin/bash
#SBATCH -A ylab
#SBATCH -J cirilong
#SBATCH -p compute
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -o /home/hxy/log/cirilong-%A_%a.out
#SBATCH -e /home/hxy/log/cirilong-%A_%a.out
#SBATCH --mem 100G
#SBATCH --array 0-6:2

workdir=/media/work/hxy/PcircNano2021/3_cirilong/2_collapse
outdir=/media/share/node11/disk3/ylab/hxy/PM1_circ_illumima_202101
circlst=$outdir/script/circFA.lst
logdir=$workdir/log
prefix=PM1

genome=/media/share/node11/disk3/ylab/hxy/reference/Talaromyces_marneffei/genome/pm1.soft.fasta
gtf=/media/share/node11/disk3/ylab/hxy/reference/Talaromyces_marneffei/annotation/pm1_scallop_final.gtf

mkdir -p $logdir

thread=$SLURM_CPUS_PER_TASK

# conda env: CIRI-long
CIRI-long collapse -i $circlst \
	-o $workdir \
	-r $gonome \
	-p $prefix \
	-a $gtf \
	-t $thread &> $logdir/${prefix}.collapes.log

# data transfor
mkdir -p $outdir/3_cirilong/2_collapse/log

cp $workdir $outdir/3_cirilong/2_collapse/
cp $logdir/${prefix}.collapes.log $outdir/3_cirilong/2_collapse/log/

# work space release
#rm -r $workdir
#rm $logdir/${prefix}.collapes.log
