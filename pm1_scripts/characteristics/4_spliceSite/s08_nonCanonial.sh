#!/usr/bin/bash

#gtf=../2_annotation/IN_Quant.clean.gtf
#
#cat all.motifRev.tsv | grep "AGAT" | cut -f5 | tr '[:upper:]' '[:lower:]' > non_canonical/aggt.tsv
#cat all.motifRev.tsv | grep "ACGT" | cut -f5 | tr '[:upper:]' '[:lower:]' > non_canonical/acga.tsv
#
#
#rm non_canonical/aggt.gtf non_canonical/acga.gtf
#
#cat non_canonical/aggt.tsv | while read id
#do
#	cat $gtf | grep $id >> non_canonical/aggt.gtf
#done
#
#cat non_canonical/acga.tsv | while read id
#do
#        cat $gtf | grep $id >> non_canonical/acga.gtf
#done

cat non_canonical/aggt.gtf | awk '{print $14}' | sed 's/;//g;s/"//g;s/\.00//g' | sort | uniq | tr 'A' 'T' > non_canonical/aggt.gene
cat non_canonical/acga.gtf | awk '{print $14}' | sed 's/;//g;s/"//g;s/\.00//g' | sort | uniq | tr 'A' 'T' > non_canonical/acga.gene

cat non_canonical/aggt.gene | while read id
do
	cat pm1_exon/pm1_intron.motifRev.tsv | grep $id >> non_canonical/aggt.motif
done

cat non_canonical/acga.gene | while read id
do
        cat pm1_exon/pm1_intron.motifRev.tsv | grep $id >> non_canonical/acga.motif
done
