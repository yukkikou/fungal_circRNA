#!/usr/bin/Rscript

library(tidyverse)

# circ info
message("*** circ_table ***")
circ_table = read_tsv("../../2_annotation/circ_gene.hyphen.map", col_names = c("circ_id","gene_id","circtrans_id"))
head(circ_table)

# circ length
message("*** circ_length ***")
circ_length = read_tsv("../../2_annotation/circ_IN.hyphen.length", col_names = c("circ_id", "ori_type", "len"))
head(circ_length)

# circ marker
message("*** circ_marker_table ***")
circ_marker_dir = "../4_clasf"
circ_marker_list = list.files(circ_marker_dir, pattern = "*grp*")
circ_marker_table = data.frame(matrix(ncol = 4, nrow = 0)) %>%
    `colnames<-`(c("circtrans_id", "transcript_id", "compose", "circ_marker"))
for (i in 1:length(circ_marker_list)){
    var = str_split(str_split(circ_marker_list[i], pattern = "-")[[1]][2], pattern = ".circ")[[1]][1]
    marker = str_split(circ_marker_list[i], pattern = "-")[[1]][1]
    
    f_tmp = read_tsv(paste(circ_marker_dir, circ_marker_list[i], sep = "/"), col_names = c("circtrans_id", "transcript_id", "compose")) %>%
        mutate(circ_marker = marker)
    head(f_tmp)
    assign(var, f_tmp)
    circ_marker_table = rbind(circ_marker_table, f_tmp)
}

head(circ_marker_table)
write_tsv(circ_marker_table, "../4_clasf/circ_marker_table.tsv")

# circ splice site
message("*** circ_ssite ***")
circ_ssite = read_tsv("../../4_spliceSite/all.motifRev.hyphen.tsv") %>%
    rename(circ_id = CIRCID,
            motif = MOTIFREV,
            strand = STRAND) %>%
    select(circ_id, motif, strand)
head(circ_ssite)

# circ flank
message("*** circ_flk ***")
circ_flkdown_motif = read_tsv("../../6_flanking/3_motif/1_streme/flk.down/flk.down.circ.hyphen.motif", col_names = c("flk_down_motif", "circ_id"))
circ_flkup_motif = read_tsv("../../6_flanking/3_motif/1_streme/flk.up/flk.up.circ.hyphen.motif", col_names = c("flk_up_motif", "circ_id"))
head(circ_flkup_motif)
head(circ_flkdown_motif)


# circ expression
message("*** circ_cpm ***")
circ_cpm = read_tsv("../../7_expression/1_circ/IN_merge_circRNA_cpm.hyphen.tsv") %>%
    pivot_longer(cols = IN_M1:IN_Y4,
                names_to = "sample",
                values_to = "CPM")
head(circ_cpm)

# merge to one table
message("*** start merge... ***")
circ_table_com = circ_table %>%
    full_join(circ_marker_table, by = c("circ_id" = "circtrans_id")) %>%
    left_join(circ_length, by = c("circ_id" = "circ_id")) %>%
    left_join(circ_ssite, by = c("circ_id" = "circ_id")) %>%
    left_join(circ_flkup_motif, by = c("circ_id" = "circ_id")) %>%
    left_join(circ_flkdown_motif, by = c("circ_id" = "circ_id")) %>%
    mutate(duplicate = duplicated(circ_id))

circ_table_cpm = circ_table %>%
    full_join(circ_marker_table, by = c("circ_id" = "circtrans_id")) %>%
    left_join(circ_length, by = c("circ_id" = "circ_id")) %>%
    left_join(circ_ssite, by = c("circ_id" = "circ_id")) %>%
    left_join(circ_flkup_motif, by = c("circ_id" = "circ_id")) %>%
    left_join(circ_flkdown_motif, by = c("circ_id" = "circ_id")) %>%
    mutate(duplicate = duplicated(circ_id)) %>%
    left_join(circ_cpm, by = c("circ_id" = "circ_id"))

head(circ_table)
write_tsv(circ_table_com, "circ_allinone.tsv")
write_tsv(circ_table_cpm, "circ_allinone_cpm.tsv")

