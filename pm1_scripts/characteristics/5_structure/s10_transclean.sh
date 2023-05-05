#!/usr/bin/bash

# uniq to one transcript

cat circ_structure.hyphen.tsv | cut -f1 | sort | uniq > circ_strcuture_all.hyphen.id
#cat circ_structure.tsv | grep "\\.01" | cut -f1 | sort | uniq | wc -l
cat circ_structure.hyphen.tsv | grep "\\.01" | cut -f1 | sort | uniq > circ_structure_first.hyphen.id
cat circ_structure.hyphen.tsv | grep "\\.01" > circ_structure_first.hyphen.tsv

#comm -13 circ_structure_first.id circ_strcuture_all.id
#
#echo `les circ_structure.tsv | grep "TM080055.c1" | head -n 1` >> circ_structure_nonfirst.tsv
#echo `les circ_structure.tsv | grep "TM020835.c1" | head -n 1` >> circ_structure_nonfirst.tsv
#echo `les circ_structure.tsv | grep "TM020255.c1" | head -n 1` >> circ_structure_nonfirst.tsv

comm -13 circ_structure_first.hyphen.id circ_strcuture_all.hyphen.id

#pm1_02:2290067-2290525
#pm1_02:716684-717101
#pm1_08:160267-160748

echo `cat circ_structure.tsv | grep "pm1_02:2290067-2290525" | head -n 1` >> circ_structure_nonfirst.hyphen.tsv
echo `cat circ_structure.tsv | grep "pm1_02:716684-717101" | head -n 1` >> circ_structure_nonfirst.hyphen.tsv
echo `cat circ_structure.tsv | grep "pm1_08:160267-160748" | head -n 1` >> circ_structure_nonfirst.hyphen.tsv

cat circ_structure_first.hyphen.tsv circ_structure_nonfirst.hyphen.tsv > circ_structure_uniq.hyphen.tsv
