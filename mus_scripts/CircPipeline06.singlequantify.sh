#!/bin/bash
#SBATCH -A ylab
#SBATCH -J singQuan
#SBATCH -p single-em
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.err
#SBATCH --mem 120G
#SBATCH --array 0-2:2

yml=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/script/mus.yml
samlst=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/0_data/CRA003317/illu_sample.list
cldir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/2_dedup
circ=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/6_cirilong/2_collapse/CIRIlongforCIRIquant.txt

workdir=/media/work/hxy/CIRIlong_zfq/7_quantify/2_singleNano
outdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/7_quantify/2_singleNano

thread=$SLURM_CPUS_PER_TASK

quandir=$workdir
logdir=$quandir/log

mkdir -p $logdir

i=$SLURM_ARRAY_TASK_ID
circ_array=(`cat $samlst`)
sample_array=(${circ_array[@]:$i:2})
clean_prefix=${sample_array[1]}

cl1=$cldir/${clean_prefix}/${clean_prefix}_R1.paired.dedup.fq.gz
cl2=$cldir/${clean_prefix}/${clean_prefix}_R2.paired.dedup.fq.gz

#conda activate CIRIquant
conda info --env
#CIRIquant
quan=$quandir/${clean_prefix}
mkdir -p $quan

echo "CIRIquant START: `date`" && CIRIquant -t $thread -v \
	-1 $cl1 \
	-2 $cl2 \
	--config $yml \
	-l 2 \
	-p ${clean_prefix}_Quant \
	-o $quan \
	--bed $circ &> $logdir/${clean_prefix}.ciriquant.log

# data transfor
mkdir -p $outdir/${clean_prefix}/circ
mkdir -p $outdir/log

out=$outdir/${clean_prefix}

cp $quan/${clean_prefix}_Quant.* $out
cp $quan/circ/${clean_prefix}_Quant_index.fa $quan/circ/${clean_prefix}_Quant.ciri $out/circ/
cp -r $quan/gene $out/
cp $logdir/${clean_prefix}.ciriquant.log $outdir/log

# work space release
rm -r $quan
