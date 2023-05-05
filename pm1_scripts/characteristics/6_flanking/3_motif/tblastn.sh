#!/bin/bash
#SBATCH -A ylab 
#SBATCH -J tblastn
#SBATCH -p compute
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.err
#SBATCH --mem 80G

pm1db=/media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/db/PM1
workdir=/media/work/hxy/blast

spdir=/media/share/node10/disk3/ylab/hxy/PM1_circ_nanopore_202101/5_BSJsite/1_quantify/3_All_l2/stat/6_flanking/3_motif/3_homologyGene/0_seq
outdir=/media/share/node10/disk3/ylab/hxy/PM1_circ_nanopore_202101/5_BSJsite/1_quantify/3_All_l2/stat/6_flanking/3_motif/3_homologyGene/1_result

thread=$SLURM_CPU_PER_TASK

# body
mkdir -p $workdir
export SINGULARITY_BIND="/media:/media,/home/hxy:/home/hxy"

# prot
singularity exec ~/singularity/blast-2.13.0/blast-2.13.0.sif tblastn -query $spdir/gene.fasta -db $pm1db \
        -evalue 1e-3 \
        -num_alignments 20 \
        -num_threads 30 \
        -out $workdir/homoGene.rlt \
        -outfmt "6 qseqid qlen qstart qend sseqid slen sstart send evalue length pident qseq sseq" \
        &> $workdir/tn.log 

# transfor
mkdir -p $outdir

cp $workdir/homoGene.rlt $workdir/tn.log $outdir

# release work space
rm $workdir/homoGene.rlt $workdir/tn.log
