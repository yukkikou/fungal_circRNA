#!/usr/bin/bash

circ=../3_bed13/circ_structure_uniq.hyphen.onetrans.tsv 
multi=multi.circ

#### single-exon gene, 51
# grp1: A, 15
cat $circ | grep "exon_0E" | awk '/o_|_o/ {print $0}' > single_exon_over-grp1_A_16.circ
# grp2: B, 36
cat $circ | grep "exon_0E" | awk '/i_.*_i/ {print $0}' > single_exon_internal-grp2_B_36.circ

#### multi-exon gene, 1830

cat $circ | grep -v "exon_0E" > $multi
cat $multi | grep 'o-o' > multi-multi.circ
cat $multi | grep "TM050088.c1" >> multi-multi.circ
cat $multi | grep -v 'o-o' |  grep -v "TM050088.c1" > multi-single.circ

### single feature
## intron
cat $multi | grep -v 'o-o' | grep "intron"  | wc -l
# intron-match
# grp3: C, 19
cat $multi | grep -v 'o-o' | awk '/m_|_m/ {print $0}' | grep "intron" > multi_single_intron_match-grp3_C_19.circ
# intron-over
cat $multi | grep -v 'o-o' | awk '/o_|_o/ {print $0}' | grep "intron" > multi_single_intron_over.circ
# intron-internal
# grp4: D, 3
cat $multi | grep -v 'o-o' | awk '/i_.*_i/ {print $0}' | grep "intron" > multi_single_intron_internal-grp4_D_3.circ


## exon
## TM050088.c1
#cat $multi | grep -v 'o-o' | grep "exon" | grep -v "TM050088.c1" | wc -l
#cat $multi | grep -v 'o-o' | grep "exon" | grep -v "TM050088.c1" | grep "[0-9]E" | wc -l
#cat $multi | grep -v 'o-o' | grep "exon" | grep -v "TM050088.c1" | grep -v "[0-9]E" | wc -l

# end-exon-match
#cat $multi | grep -v 'o-o' | grep "exon" | grep -v "TM050088.c1" | grep "[0-9]E" | awk '/m_|_m/ {print $0}' > multi_single_end-exon_match.circ
# end-exon-over
#cat $multi | grep -v 'o-o' | grep "exon" | grep -v "TM050088.c1" | grep "[0-9]E" | awk '!/m_|_m/ {print $0}' | awk '/o_|_o/ {print $0}' > multi_single_end-exon_over.circ
# end-exon-internal
#cat $multi | grep -v 'o-o' | grep "exon" | grep -v "TM050088.c1" | grep "[0-9]E" | awk '!/m_|_m/ {print $0}' | awk '!/o_|_o/ {print $0}' > multi_single_end-exon_internal.circ

# internal-exon exact match
# grp5: E, 294
cat $multi | grep -v 'o-o' | grep "exon" | grep -v "TM050088.c1" | grep -v "[0-9]E" | awk '/m_.*_m/ {print $0}' > multi_single_internal_exon_exact_match-grp5_E_296.circ
# grp6: F, 694
cat $multi | grep -v 'o-o' | grep "exon" | grep -v "TM050088.c1" | awk '!/m_.*_m/ {print $0}' > multi_single_exons-grp6_F_702.circ

### multi feature
#cat $multi | grep 'o_o' | wc -l 

# intron edge 

# grp7: G, 42
cat multi-multi.circ | awk '!/.*_intron_[0-9]E?_[imo]$/ {print $0}' | awk '/\t[iom]_intron_*/ {print $0}' > multi_intron_start-grp7_G_40.circ
# grp8: H, 57
cat multi-multi.circ | awk '!/\t[iom]_intron_*/ {print $0}' | awk '/.*_intron_[0-9]E?_[imo]$/ {print $0}' > multi_intron_end-grp8_H_57.circ
# grp9: I, 2
cat multi-multi.circ | awk '/.*_intron_[0-9]E?_[imo]$/ {print $0}' | awk '/\t[iom]_intron_*/ {print $0}' > multi_intron_start_end-grp9_2.circ
# grep10: J, 721
cat multi-multi.circ | awk '!/\t[iom]_intron_*/ {print $0}' | awk '!/.*_intron_[0-9]E?_[imo]$/ {print $0}' > multi_exon_exon-grp10_725.circ
# special
# grp10: J1, 7
for i in {i-i,i-o,i-m,m-m,m-i,m-o,o-i,o-m}
do
    echo "*** $i ***"
    cat $circ | grep $i | sort | uniq
    cat $circ | grep $i | sort | uniq >> multi_exon_eejunc-grp10_J1_6.circ
done
