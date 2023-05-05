#!/bin/bash

mkdir -p stat && cd stat

for sm in {C1,C2,EC1,EC3,Y2,Y3,EY3,EY4,All}
do
	cat ${sm}.motifRev.txt | cut -f10 | sed '1d' | tr [:lower:] [:upper:] | sort | uniq -c | tr -s " " > ${sm}.stat.txt
done
