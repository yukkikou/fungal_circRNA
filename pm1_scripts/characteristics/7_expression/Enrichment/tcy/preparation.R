omicsbox_table <- read.delim("blast2go/blast2go_annot.annot",
                             row.names=NULL,header = F,
                             col.names = c('GENEID','ANNOTATION','DESCRIPTION'),
                             quote="",
                             stringsAsFactors=FALSE)

PM1go2term<-go2term(omicsbox_table$ANNOTATION)
rownames(PM1go2term)<-PM1go2term$go_id
PM1go2ont<-go2ont(omicsbox_table$ANNOTATION)
rownames(PM1go2ont)<-PM1go2ont$go_id
PM1GO<-data.frame(row.names =PM1go2term$go_id,
                  go_id=PM1go2term$go_id,
                  term=paste0(PM1go2term$Term,' (',PM1go2term$go_id,')'),
                  ont=PM1go2ont[PM1go2term$go_id,2],
                  stringsAsFactors = F)

omicsbox_table<-omicsbox_table[omicsbox_table$ANNOTATION%in%PM1GO$go_id,]
PM1GO<-data.frame(Gene=paste0(substr(omicsbox_table$GENEID,1,8),'.00'),
                  Transcript=omicsbox_table$GENEID,
                  GO=omicsbox_table$ANNOTATION,
                  TERM=PM1GO[omicsbox_table$ANNOTATION,'term'],
                  Ontology=PM1GO[omicsbox_table$ANNOTATION,'ont'],
                  # GOlevel=omicsbox_table$Level,
                  stringsAsFactors = F)


write.table(PM1GO,'pm1.go',col.names = T,quote = F,
            row.names = F,sep = '\t')
PM1GOgene<-PM1GO[,-2]
PM1GOgene<-PM1GOgene[!duplicated(PM1GOgene),]

PM1db<-sapply(rownames(dat_de),
              function(x){
                return(PM1GOgene$GO[PM1GOgene$Gene==x])})


GOlevel <- read.delim("blast2go/blast2go_go_propagation.txt",
                      row.names=NULL,header = T,
                      # col.names = c('GENEID','ANNOTATION','DESCRIPTION'),
                      quote="",
                      stringsAsFactors=FALSE)
GOlevel<-GOlevel[,c(2,4)]
GOlevel<-GOlevel[!duplicated(GOlevel),]
GOlevel<-GOlevel[GOlevel$Level>3,]
GOlevel<-GOlevel[GOlevel$Level<6,]
GOlevel<-unique(GOlevel$GO)


myCollection<-GOCollection(
  as.character(
    unique(
      PM1GOgene$GO)))
slim <- getOBOCollection('blast2go/goslim_yeast.obo')
slim.go<-slim@ids

rownames(GOlevel)<-GOlevel$GO


PM1db<-sapply(rownames(dds),
              function(x){
                return(PM1GOgene$GO[PM1GOgene$Gene==x])})