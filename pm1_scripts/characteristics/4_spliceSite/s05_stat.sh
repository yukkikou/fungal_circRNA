#!/usr/bin/bash

tail -n+2 all.motifRev.txt | awk '{print $10}' | tr '[:lower:]' '[:upper:]' | sort | uniq -c | tr -s " " | awk -v OFS='\t' '{print $2,$1,$1/3567}' | sort -k2 -n > all.motifRev.table

awk '{sum+=$3} END {print sum}' all.motifRev.table
