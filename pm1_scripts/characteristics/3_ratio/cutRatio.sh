#!/usr/bin/bash

for i in `ls 0_separate`
do
	tail -n+6 0_separate/$i | awk '{print $18}' | sed 's/;//g' > $(basename $i ".gtf").ratio
done
