#!/bin/bash
#SBATCH -A ylab
#SBATCH -J circZFQ
#SBATCH -p compute
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.err
#SBATCH --mem 80G
#SBATCH --array 0-6:2

fqdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/0_data/CRA003317
samlst=$fqdir/illu_sample.list

workdir=/media/work/hxy/CIRIlong_zfq/1_clean
outdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/1_clean

i=$SLURM_ARRAY_TASK_ID
thread=$SLURM_CPUS_PER_TASK

logdir=$workdir/log
repdir=$workdir/report

mkdir -p $logdir $repdir

circ_array=(`cat $samlst`)
sample_array=(${circ_array[@]:$i:2})

raw_prefix=${sample_array[0]}
clean_prefix=${sample_array[1]}

fq1=$fqdir/${raw_prefix}/${raw_prefix}_f1.fq.gz
fq2=$fqdir/${raw_prefix}/${raw_prefix}_r2.fq.gz
cl1=$workdir/${clean_prefix}_R1.paired.fq.gz
cl2=$workdir/${clean_prefix}_R2.paired.fq.gz

# conda env: QC
fastp --in1 $fq1 --in2 $fq2 \
	--out1 $cl1 --out2 $cl2 \
	--detect_adapter_for_pe \
        --cut_tail --cut_front --cut_mean_quality 20 \
        -q 20 --length_required 140 --length_limit 160 \
        -w $thread \
        -h $repdir/${clean_prefix}_fastp.html \
        -j $repdir/${clean_prefix}_fastp.json &> $logdir/${clean_prefix}.fastp.log

# data transfor
mkdir -p $outdir/{log,report}
mkdir -p $outdir/${clean_prefix}

cp $cl1 $cl2 $outdir/${clean_prefix}
cp $repdir/${clean_prefix}_fastp.* $outdir/report/
cp $logdir/${clean_prefix}.fastp.log $outdir/log/

# work space release
rm ${clean_prefix}_* $repdir/${clean_prefix}_fastp.* $logdir/${clean_prefix}.fastp.log
