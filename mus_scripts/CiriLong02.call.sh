#!/bin/bash
#SBATCH -A ylab
#SBATCH -J cl_ZFQ
#SBATCH -p single-em
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.err
#SBATCH --mem 100G
#SBATCH --array 0-2:2

# mamba env: CIRI-long
# 30 threads took over a whole day for each task
genome=/media/share/node10/disk3/ylab/hxy/reference/Mus_musculus/Genome/GRCm38.p4.genome.fa
gtf=/media/share/node10/disk3/ylab/hxy/reference/Mus_musculus/Annotation/gencode.vM10.annotation.gtf

workdir=/media/work/hxy/CIRIlong_zfq/6_cirilong/1_call
cldir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/5_nano
outdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/6_cirilong/1_call
samlst=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/0_data/CRA003317/nano_sample.list

mkdir -p $workdir

i=$SLURM_ARRAY_TASK_ID
thread=$SLURM_CPUS_PER_TASK

circ_array=(`cat $samlst`)
sample_array=(${circ_array[@]:$i:2})
clean_prefix=${sample_array[1]}
cl=$cldir/${clean_prefix}.clean.fq

out=$workdir/${clean_prefix}
logdir=$workdir/log

mkdir -p $out $logdir

# Index
#cd $genomedir
#bwa index -a bwtsw $genome

# conda env: CIRI-long
CIRI-long call -i $cl \
	-o $out \
	-r $genome \
	-p ${clean_prefix} \
	-a $gtf \
	-t $thread &> $logdir/${clean_prefix}.call.log

# data transfor
mkdir -p $outdir/log

cp -r $out $outdir
cp $logdir/${clean_prefix}.call.log $outdir/log

# work space release
rm -r $out
rm $logdir/${clean_prefix}.call.log
