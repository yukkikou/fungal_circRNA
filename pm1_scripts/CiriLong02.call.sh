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

# 30 threads took over a whole day for each task
workdir=/media/work/hxy/PcircNano2021/3_cirilong/1_call
outdir=/media/share/node11/disk3/ylab/hxy/PM1_circ_illumima_202101
samlst=$outdir/script/circ_sample.list

genome=/media/share/node11/disk3/ylab/hxy/reference/Talaromyces_marneffei/genome/pm1.soft.fasta
gtf=/media/share/node11/disk3/ylab/hxy/reference/Talaromyces_marneffei/annotation/pm1_scallop_final.gtf

mkdir -p $workdir

i=$SLURM_ARRAY_TASK_ID
thread=$SLURM_CPUS_PER_TASK

circ_array=(`cat $samlst`)
sample_array=(${circ_array[@]:$i:2})
clean_prefix=${sample_array[1]}
cl=$outdir/1_clean/${clean_prefix}.clean.fq

out=$workdir/${clean_prefix}
logdir=$workdir/log

mkdir -p $out $logdir

# conda env: CIRI-long
CIRI-long call -i $cl \
	-o $out \
	-r $gonome \
	-p ${clean_prefix} \
	-a $gtf \
	-t $thread &> $logdir/${clean_prefix}.call.log

# data transfor
mkdir -p $outdir/3_cirilong/1_call/log

cp $out $outdir/3_cirilong/1_call/
cp $logdir/${clean_prefix}.call.log $outdir/3_cirilong/1_call/log/

# work space release
#rm -r $out
#rm $logdir/${clean_prefix}.call.log
