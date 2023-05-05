cat circlong_suppIN.outv2 | sed 's/:/\t/g' | awk '{if($2!=$7&&$3!=1) print $0}' | awk -v OFS='' '{print $1,":",$2}' > unAnnotation.id
cat circlong_suppIN.outv2 | sed 's/:/\t/g' | awk '{if($2!=$7&&$3==1) print $0}' | awk -v OFS='' '{print $1,":",$2}' > unAnnoSingle.id
cp circlong_suppIN.bed12 circlong_suppIN.bed12.bak
cat circlong_suppIN.bed12 | sed 's/,/_/g' > circlong_suppIN.r.bed12
