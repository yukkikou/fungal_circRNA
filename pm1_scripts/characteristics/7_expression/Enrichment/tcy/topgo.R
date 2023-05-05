#TopGO####
hvg<-list('dem2y'=unlist(sapply(1:10,function(x){
  df<-res[,strsplit2(colnames(dat_de),'_vs_')[,1]==paste0('Time',x)|
            strsplit2(colnames(dat_de),'_vs_')[,2]==paste0('Time',x)]
  if(x>10){
    idx<-c(rep(1,x-11),rep(-1,21-x))
    }else{idx<-c(rep(1,x-1),rep(-1,11-x))}
  rownames(dat_de)[apply(df*idx>1,1,sum)>8|apply(-1>df*idx,1,sum)>8]
  })),
  'dey2m'=unlist(sapply(11:20,function(x){
    df<-res[,strsplit2(colnames(dat_de),'_vs_')[,1]==paste0('Time',x)|
              strsplit2(colnames(dat_de),'_vs_')[,2]==paste0('Time',x)]
    if(x>10){
      idx<-c(rep(1,x-11),rep(-1,21-x))
    }else{idx<-c(rep(1,x-1),rep(-1,11-x))}
    rownames(dat_de)[apply(df*idx>1,1,sum)>8|apply(-1>df*idx,1,sum)>8]
  })),
  'tao_m2y'=tao_m2y,
  'tao_y2m'=tao_y2m,
  'HighVariableGene_m2y'=HighVariableGene_m2y,
  'HighVariableGene_y2m'=HighVariableGene_y2m)

goresult<-as.data.frame(
  t(
    data.frame(row.names = c("GO.ID","Term","Annotated","Significant",
                             "Expected","classic",'classic.p.adjust',
                             "weightFisher",'weightFisher.p.adjust',
                             'Ontology','type'))
  )
)


for (cl in 1:6) {
  
  for (ont in c('BP','MF','CC')) {
    # dat_de<-dat_de_total
    genelist<-factor(as.integer(rownames(dat_de)%in%
                                  hvg[[cl]]))
    names(genelist)<-rownames(dat_de)
    
    GOdata <- new("topGOdata",
                  ontology = ont,
                  allGenes = genelist,
                  annot = annFUN.gene2GO,
                  nodeSize = 1,
                  gene2GO = PM1db)
    resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    weight01.fisher <- runTest(GOdata, statistic = "fisher")
    allRes <- GenTable(GOdata,
                       classic = resultFis,
                       weightFisher = weight01.fisher,
                       orderBy = "classic",
                       ranksOf = "classic",
                       topNodes = length(usedGO(GOdata)))
    go<-data.frame(allRes,stringsAsFactors = F)
    go$classic<-as.numeric(gsub('<','',go$classic))
    go$weightFisher<-as.numeric(gsub('<','',go$weightFisher))
    go$classic.p.adjust<-p.adjust(go$classic)
    go$weightFisher.p.adjust<-p.adjust(go$weightFisher)
    go$Ontology<-rep(ont,dim(go)[1])
    go$Module<-rep(cl,dim(go)[1])
    goresult<-rbind(go,goresult)
  }
}

topgo_result<-goresult

df<-topgo_result

df$Description<-paste0(df$Term,' (',
                       df$GO.ID,')')
df<-df[df$weightFisher.p.adjust<0.05,]
df<-df[df$classic.p.adjust<0.1,]
