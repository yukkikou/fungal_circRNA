#!/usr/bin/bash

gtf=IN_Quant.singlegene.clean.sorted.gtf

cat $gtf | \
	awk 'BEGIN{OFS="";oldtr="";}{
    if($14==""){
        print $0," circtrans_id ","_",$10,"_.c",0,"\";";
    } else {
		if($3=="circRNA") {
			tr=$14;
    		if(oldtr!=tr) {
	    		start=1; oldtr=tr;
		   		print $0," circtrans_id ","_",tr,"_.c",start,"\";";
		    	start+=1;
		   	} else {
		     	print $0," circtrans_id ","_",tr,"_.c",start,"\";";
		    	start+=1;
		   	} 
        } else {
			print $0;
		}
    }
}' | sed 's/_"/"/g;s/";_//g;s/\.00\./\./g' > IN_Quant.singlegene.clean.circid.sorted.gtf

