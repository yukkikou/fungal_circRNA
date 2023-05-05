#!/bin/bash
#SBATCH -A ylab
#SBATCH -J CIRI2_ZFQ
#SBATCH -p compute
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.out
#SBATCH --mem 240G
#SBATCH --array 0-6:2

# mamba env: bwa
gtf=/media/share/node10/disk3/ylab/hxy/reference/Mus_musculus/Annotation/gencode.vM10.annotation.gtf
genome=/media/share/node10/disk3/ylab/hxy/reference/Mus_musculus/Genome/GRCm38.p4.genome.fa
bwaindex=/media/share/node10/disk3/ylab/hxy/reference/Mus_musculus/Genome/GRCm38.p4.genome.fa

samlst=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/0_data/CRA003317/illu_sample.list

cldir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/2_dedup
workdir=/media/work/hxy/CIRIlong_zfq
outdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq

thread=$SLURM_CPUS_PER_TASK

bamdir=$workdir/3_bwabam
circdir=$workdir/4_ciri2
bamlogdir=$bamdir/log
circlogdir=$circdir/log

mkdir -p $bamlogdir $circlogdir

i=$SLURM_ARRAY_TASK_ID
circ_array=(`cat $samlst`)
sample_array=(${circ_array[@]:$i:2})
clean_prefix=${sample_array[1]}

cl1=$cldir/${clean_prefix}/${clean_prefix}_R1.paired.dedup.fq.gz
cl2=$cldir/${clean_prefix}/${clean_prefix}_R2.paired.dedup.fq.gz

bam=$bamdir/${clean_prefix}.sam
circ=$circdir/${clean_prefix}.circ.txt

# CIRI2
bwa mem -T 19 -t $thread -o $bam $bwaindex $cl1 $cl2 &> $bamlogdir/${clean_prefix}.bwa.log

perl /home/hxy/software/CIRI_v2.0.6/CIRI2.pl -I $bam \
	-O $circ \
	-F $genome \
	-A $gtf \
	-G $circlogdir/${clean_prefix}.ciri2.log \
	-high -T $thread

# data transfor
mkdir -p $outdir/3_bwabam/${clean_prefix}
mkdir -p $outdir/3_bwabam/log
mkdir -p $outdir/4_ciri2/${clean_prefix}
mkdir -p $outdir/4_ciri2/log

cp $bam $outdir/3_bwabam/${clean_prefix}
cp $circ $outdir/4_ciri2/${clean_prefix}
cp $bamlogdir/${clean_prefix}.bwa.log $outdir/3_bwabam/log/
cp $circlogdir/${clean_prefix}.ciri2.log $outdir/4_ciri2/log/

# work space release
rm $bam $circ $bamlogdir/${clean_prefix}.bwa.log $circlogdir/${clean_prefix}.ciri2.log
