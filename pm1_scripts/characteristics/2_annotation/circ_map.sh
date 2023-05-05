#!/usr/bin/bash

gtf=IN_Quant.singlegene.clean.circid.sorted.gtf

cat $gtf | awk -v OFS='\t' '{
    if ($16==""){
        print $10,$10,$14
    } else {
        print $10,$14,$16
    }
}' | sed 's/"//g;s/;//g' > circ_gene.sorted.vertical.map

cat $gtf | awk -v OFS='\t' '{
    if ($16==""){
        print $10,$10,$14
    } else {
        print $10,$14,$16
    }
}' | sed 's/"//g;s/;//g;s/|/-/g' > circ_gene.sorted.hyphen.map
