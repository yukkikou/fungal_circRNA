#!/usr/bin/bash

# CIRI2
ciri2dir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/4_ciri2
ciriall=IN_all_total.circ.txt
ciribed=CIRI2forCIRIquant.txt

cd $ciri2dir
cat `find -name "IN_total_rep[1-2].circ.txt"` | sort | uniq > $ciriall
cat $ciriall | grep -v "circRNA_ID" | awk -v OFS='\t' '{print $2,$3,$4,$1,".",$11}' |  sed 's/"//g;s/;//g' > $ciribed

# CIRI-long
longdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/6_cirilong/2_collapse
long=Mus.info
longbed=CIRIlongforCIRIquant.txt

cd $longdir
cat $long | awk -v OFS='\t' '{print $1,$4,$5,$10,$8,$7}' | sed 's/"//g;s/;//g;s/-/|/;s/None/n\/a/' | sort | uniq > $longbed

outdir=/media/share/node10/disk3/ylab/hxy/CIRIlong_zfq/7_quantify
outbed=AlltoCIRIquant.txt

mkdir -p $outdir && cd $outdir
cat $ciri2dir/$ciribed $longdir/$longbed | sort | uniq > $outbed
