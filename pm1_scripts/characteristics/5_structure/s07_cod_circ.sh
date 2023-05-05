#!/usr/bin/bash

gtf=circlong_suppIN.final.gtf
#sed -i 's/|/"/g;s/_/;/' $gtf

cat $gtf | uniq | \
	awk 'BEGIN{OFS="";oldtr="";}{
		if($3=="circExon") {
			tr=$14;
    		if(oldtr!=tr) {
	    		start=1; oldtr=tr;
		   		print $0," circexon_num ",start,";"," circexon_size ",$5-$4+1,";"," circexon_id ","_",tr,"_.e",start,"\";";
		    	start+=1;
		   	} else {
		     	print $0," circexon_num ",start,";"," circexon_size ",$5-$4+1,";"," circexon_id ","_",tr,"_.e",start,"\";";
		    	start+=1;
		   	} 
        } else {
			print $0;
		}
}' | sed 's/";_//g;s/_"/"/g;s/pm1;/pm1_/g' > circlong_suppIN.final.circid.gtf

isige ~/singularity/bedtools/bedtools-2.30.0.sif bedtools intersect -s -a circlong_suppIN.final.circid.gtf -b ../1_geneClas/pm1_scallop_final_num.bed -wao > circ_exon_feat.inter
cat circ_feat_exon.inter | awk -v OFS='\t' '$21!="."{print $4,$5,$14,$18,$20,$22,$23,$24,$25,$27}' | sed 's/"//g;s/;//g' > circ_feat_exon.tsv
cat circ_feat_exon.inter | awk -v  OFS='\t' '$21!="."{print $4,$5,$10,$18,$20,$22,$23,$24,$25,$27}' | sed 's/"//g;s/;//g' > circ_feat_exon.hyphen.tsv
