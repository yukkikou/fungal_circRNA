for i in `find -name "*_Quant.gtf"`
do
	tail -n +6 $i | awk '{print $10}' | sed 's/"//g;s/;//g'| sort > stat/$(basename $i ".gtf").id
done
