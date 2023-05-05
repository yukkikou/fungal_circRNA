#!/bin/bash
#SBATCH -A ylab
#SBATCH -J chop_ZFQ
#SBATCH -p single-em
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.out
#SBATCH --mem 100G
#SBATCH --array 0-2:2

# mamba env: CIRI-long
fqdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/0_data/CRA003317
genome=/media/share/node10/disk3/ylab/hxy/reference/Mus_musculus/Genome/GRCm38.p4.genome.fa

workdir=/media/work/hxy/CIRIlong_zfq
outdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq
samlst=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/0_data/CRA003317/nano_sample.list

mkdir -p $workdir/5_nano

i=$SLURM_ARRAY_TASK_ID
thread=$SLURM_CPUS_PER_TASK

circ_array=(`cat $samlst`)
sample_array=(${circ_array[@]:$i:2})

raw_prefix=${sample_array[0]}
clean_prefix=${sample_array[1]}

fq=$fqdir/${raw_prefix}/${raw_prefix}.fq.gz
cl=$workdir/5_nano/${clean_prefix}.clean.fq

# remove adapter of circfl
/home/hxy/software/Porechop/porechop-runner.py -t $thread \
	-i $fq -o $cl --barcode_threshold 95 \
	--check_reads 1000

# data transfor
mkdir -p $outdir/5_nano
cp $cl $outdir/5_nano

# work space release
rm $cl
