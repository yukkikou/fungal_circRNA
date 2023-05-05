#!/usr/bin/bash

workdir=/media/share/node10/disk3/ylab/hxy/PM1_circ_nanopore_202101/5_BSJsite/1_quantify/3_All_l2/stat/6_flanking/3_motif
flkUpMotif=$workdir/1_streme/flk.up/streme.txt
flkDownMotif=$workdir/1_streme/flk.down/streme.txt
outdir=$workdir/2_tomtom

motifDB=/media/share/node10/disk3/ylab/hxy/reference/MotifDB/motif_databases
mfile1=JASPAR/JASPAR2022_CORE_fungi_redundant_v2.meme
mfile2=RNA/Ray2013_rbp_All_Species.dna_encoded.meme
mfile3=YEAST/SwissRegulon_s_cer.meme
lst4=CISBP-RNA/_fungi_species.lst

export SINGULARITY_BIND='/media:/media,/home/hxy:/home/hxy'

#for i in {$mfile1,$mfile2,$mfile3}
#do
#    echo "*** $i ***"
#    singularity exec ~/singularity/meme-5.4.1/meme-5.4.1.sif tomtom \
#        -o $outdir/flk.up$(basename $i ".meme") \
#        $flkUpMotif $motifDB/$i &> flk.up.tomtom.log
#
#    singularity exec ~/singularity/meme-5.4.1/meme-5.4.1.sif tomtom \
#        -o $outdir/flk.down$(basename $i ".meme") \
#        $flkDownMotif $motifDB/$i &> flk.down.tomtom.log
#
#    mv flk.up.tomtom.log $outdir/flk.up$(basename $i ".meme") 
#    mv flk.down.tomtom.log $outdir/flk.down$(basename $i ".meme")
#
#    tdir=`dirname $i`
#    echo $tdir
#    mkdir -p $outdir/$tdir
#    mv $outdir/flk.up$(basename $i ".meme") $outdir/flk.down$(basename $i ".meme") $outdir/$tdir
#done


for i in `cat $motifDB/$lst4`
do
    echo "*** $i ***"
    singularity exec ~/singularity/meme-5.4.1/meme-5.4.1.sif tomtom \
        -o $outdir/flk.up$(basename $i ".meme") \
        $flkUpMotif $motifDB/$i &> flk.up.tomtom.log

    singularity exec ~/singularity/meme-5.4.1/meme-5.4.1.sif tomtom \
        -o $outdir/flk.down$(basename $i ".meme") \
        $flkDownMotif $motifDB/$i &> flk.down.tomtom.log

    mv flk.up.tomtom.log $outdir/flk.up$(basename $i ".meme") 
    mv flk.down.tomtom.log $outdir/flk.down$(basename $i ".meme")

    tdir=`dirname $i`
    tname=`basename $i ".meme"`
    echo $tdir/$tname
    mkdir -p $outdir/$tdir/$tname
    mv $outdir/flk.up$(basename $i ".meme") $outdir/flk.down$(basename $i ".meme") $outdir/$tdir/$tname
done
