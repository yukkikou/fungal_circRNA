#!/usr/bin/bash

conUpFa=exon.up.fa
conDownFa=exon.down.fa

flkUpFa=flk.up.fa
flkDownFa=flk.down.fa

export SINGULARITY_BIND='/media:/media,/home/hxy:/home/hxy'

singularity exec ~/singularity/meme-5.4.1/meme-5.4.1.sif streme --p $flkDownFa \
   --n $conDownFa \
   --order 2 \
   -o $(basename $flkDownFa ".fa") \
   --dna &> $(basename $flkDownFa ".fa").log
