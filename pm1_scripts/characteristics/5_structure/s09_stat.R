#!/usr/bin/Rscript

library(tidyverse)

test = read_tsv("circ_feat_exon.hyphen.tsv", col_names = c("ce_ss","ce_es","circtrans_id","ce_size","ce_id","e_ss0","e_es","e_id","e_size","overlap")) %>%
    mutate(e_ss = e_ss0+1)

test$e_type = str_split(test$e_id, '-', simplify = T)[,1]
test$e_trans = str_split(test$e_id, '-', simplify = T)[,2]
test$e_num = str_split(test$e_id, '-', simplify = T)[,3]

# ss_in/match/over, es_in/match/over
for (i in 1:dim(test)[1]){
    test[i, "ss_border"] = ifelse(test[i,"ce_ss"] < test[i, "e_ss"], "o",
                        ifelse(test[i, "ce_ss"] == test[i, "e_ss"], "m", "i"))
    test[i, "es_border"] = ifelse(test[i,"ce_es"] < test[i, "e_es"], "i",
                        ifelse(test[i, "ce_es"] == test[i, "e_es"], "m", "o"))
}


test_cod = test %>% mutate(
    ce_pos = paste(ss_border, e_type, e_num, es_border, sep = "_")
) %>% filter(e_type != 'CDS')

# group_by(circtrans_id, e_trans, e_type) %>%

circtrans_id = unique(test_cod$circtrans_id)
test_merge = data.frame(matrix(nrow = 0, ncol = 3)) %>%
    `colnames<-`(c('circtrans_id','e_trans','ce_cod'))
n = 1
for (i in circtrans_id){
    circ_tmp = test_cod %>% filter(circtrans_id == i)
    trans_id_tmp = sort(unique(circ_tmp$e_trans))
    for (j in trans_id_tmp){
        circ_trans_tmp = circ_tmp %>% filter(e_trans == j) %>%
            arrange(e_ss)
        ce_cod = paste(circ_trans_tmp['ce_pos'], sep = "-")
        test_merge[n,'circtrans_id'] = i
        test_merge[n,'e_trans'] = j
        test_merge[n,'ce_cod'] = ce_cod
        test_merge[n, 'ce_cod_merge'] = str_c(unlist(circ_trans_tmp['ce_pos']), collapse = "-")
        n = n + 1
    } 
}

write_tsv(test_merge[,c(1,2,4)], "circ_structure.hyphen.tsv")
