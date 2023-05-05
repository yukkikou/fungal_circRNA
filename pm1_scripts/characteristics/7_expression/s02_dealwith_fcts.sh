#!/usr/bin/bash

workdir=/media/data5/hxy/PM1_linear_illumina_2021/3_featureCounts

#bothGeneFcts=$workdir/both/gene_count_bothend.tsv
#bothTransFcts=$workdir/both/transcript_count_bothend.tsv
bothGeneFcts=$workdir/gene_count.tsv
bothTransFcts=$workdir/transcript_count.tsv


mkdir -p $workdir/results

tail -n+2 $bothGeneFcts \
	| cut -f 1,7- \
	| sed 's#/media/data5/hxy/PM1_linear_illumina_2021/2_bam/##g' \
	| sed 's#\.hisat2\.sorted\.bam##g' > $workdir/results/$(basename $bothGeneFcts ".tsv").clean.tsv

tail -n+2 $bothTransFcts \
	| cut -f 1,7- \
	| sed 's#/media/data5/hxy/PM1_linear_illumina_2021/2_bam/##g' \
	| sed 's#\.hisat2\.sorted\.bam##g' > $workdir/results/$(basename $bothTransFcts ".tsv").clean.tsv


tail -n +2 $bothGeneFcts \
    | cut -f 1,6 > $workdir/results/$(basename $bothGeneFcts ".tsv").length

tail -n +2 $bothTransFcts \
    | cut -f 1,6 > $workdir/results/$(basename $bothTransFcts ".tsv").length


