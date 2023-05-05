#!/bin/bash
#SBATCH -A ylab
#SBATCH -J porechop
#SBATCH -p compute
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -o /home/hxy/log/porechop-%A_%a.out
#SBATCH -e /home/hxy/log/porechop-%A_%a.out
#SBATCH --mem 100G
#SBATCH --array 0-6:2

fqdir=/media/share/node13/disk3/ylab/hxy/circFL_202101/1_fq/
genomedir=/media/share/node11/disk3/ylab/hxy/reference/Talaromyces_marneffei/genome
genome=pm1.soft.fasta

workdir=/media/work/hxy/PcircNano2021
outdir=/media/share/node11/disk3/ylab/hxy/PM1_circ_illumima_202101
samlst=$outdir/script/circ_sample.list

mkdir -p $workdir/1_clean

i=$SLURM_ARRAY_TASK_ID
thread=$SLURM_CPUS_PER_TASK

circ_array=(`cat $samlst`)
sample_array=(${circ_array[@]:$i:2})

raw_prefix=${sample_array[0]}
clean_prefix=${sample_array[1]}

fq=$fqdir/${raw_prefix}.pass.fastq
cl=$workdir/1_clean/${clean_prefix}.clean.fq

# bulid index
cd $genomedir
bwa index -a bwtsw $genome

# remove adapter of circfl
/home/hxy/software/Porechop/porechop-runner.py -t $thread \
	-i $fq -o $cl --barcode_threshold 95 \
	--check_reads 1000

# data transfor
mkdir -p $outdir/1_clean
cp $cl $outdir/1_clean

# work space release
rm $cl
