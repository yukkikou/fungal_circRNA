#!/bin/bash
#SBATCH -A ylab
#SBATCH -J cDE
#SBATCH -p compute
#SBATCH -D /media/work/hxy/
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -o /home/hxy/log/%x-%A_%a.out
#SBATCH -e /home/hxy/log/%x-%A_%a.err
#SBATCH --mem 120G


CIRIsp=/media/share/node10/disk3/ylab/hxy/PM1_circ_nanopore_202101/script/sampleDE.list
Stringtiesp=/media/share/node10/disk3/ylab/hxy/PM1_circ_nanopore_202101/script/samplePrepDE.list
#workdir=/media/work/hxy/PM1_circ_nanopore_202101/5_BSJsite
outdir=/media/share/node10/disk3/ylab/hxy/PM1_circ_nanopore_202101/5_BSJsite/1_quantify/3_All_l2/stat/7_expression

prefix=$outdir/IN_merge

thread=$SLURM_CPUS_PER_TASK
logdir=$outdir/log

mkdir -p $logdir

# CIRIquant
# step1
prep_CIRIquant -i $CIRIsp \
    --lib ${prefix}_library_info.csv \
    --circ ${prefix}_circRNA_info.csv \
    --bsj ${prefix}_circRNA_bsj.csv \
    --ratio ${prefix}_circRNA_ratio.csv &> $logdir/prep_CIRIquant.log

# step2
python ~/script/5_stringtie/prepDE.py -i $Stringtiesp \
    -g ${prefix}_gene_count_matrix.csv \
    -t ${prefix}_transcript_count_matrix.csv &> $logdir/prepDE.log

# step3
CIRI_DE_replicate --lib ${prefix}_library_info.csv \
    --bsj ${prefix}_circRNA_bsj.csv \
    --gene ${prefix}_gene_count_matrix.csv \
    --out ${prefix}_circRNA_de.tsv &> $logdir/CIRI_DE_replicate.log

# others expression matrix
python ~/script/5_stringtie/getTPM.py -i $Stringtiesp \
    -g ${prefix}_gene_TPM_matrix.csv \
    -t ${prefix}_transcript_TPM_matrix.csv &> $logdir/getTPM.log

python ~/script/5_stringtie/getFPKM.py -i $Stringtiesp \
    -g ${prefix}_gene_FPKM_matrix.csv \
        -t ${prefix}_transcript_FPKM_matrix.csv &> $logdir/getFPKM.log
