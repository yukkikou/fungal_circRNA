#!/usr/bin/bash

cat circlong_suppIN.r.clean.bed12 | sed 's/|/,/g' > circlong_suppIN.clean.bed12

sige ~/singularity/bedtools/bedtools-2.30.0.sif bedtools getfasta -s -name -split -fi /media/share/node10/disk3/ylab/hxy/reference/Talaromyces_marneffei/genome/pm1.soft.fasta -bed circlong_suppIN.clean.bed12 > circlong_suppIN.circ.fasta
