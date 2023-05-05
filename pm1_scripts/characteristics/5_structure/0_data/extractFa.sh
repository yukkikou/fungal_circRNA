#!/usr/bin/bash

fa=/media/share/node10/disk3/ylab/hxy/PM1_circ_nanopore_202101/3_cirilong/2_collapse/PM1.reads

idlst=circlong_suppIN_id

rm circlong_suppIN.lst

time cat $idlst | while read id
do
	cat $fa | grep $id >> circlong_suppIN.lst
done

