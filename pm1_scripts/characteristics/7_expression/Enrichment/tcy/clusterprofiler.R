require(doParallel)
require(foreach)
stopImplicitCluster()
#registerDoParallel(104)
goresult<-foreach(idx=1:max(info.norm$tclust),
                  # idx=1:12,
                  .combine = 'rbind',
                  .packages = c('clusterProfiler','rvest','httr'))%dopar%
  {
    genelist<-info.norm$gene[info.norm$tclust==idx]
    go<-as.data.frame(
      t(
        data.frame(row.names = c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust",
                                 "qvalue","geneID","Count"))
      )
    )
    go[1,]=1
    go_bp<-enricher(genelist,
                    pAdjustMethod = 'fdr',minGSSize = 0,
                    # TERM2GENE=PM1GOgene[PM1GOgene$Ontology=='BP',c(2,1)],
                    # TERM2NAME=PM1GOgene[PM1GOgene$Ontology=='BP',c(2,3)],
                    TERM2GENE=PM1GOgene[PM1GOgene$Gene%in%rownames(tclust_cts)&PM1GOgene$Ontology=='BP',c(2,1)],
                    TERM2NAME=PM1GOgene[PM1GOgene$Gene%in%rownames(tclust_cts)&PM1GOgene$Ontology=='BP',c(2,3)],
                    pvalueCutoff=1,
                    qvalueCutoff = 1)
    if(is.null(go_bp)){
      go_bp<-go
    }else{
      go_bp<-data.frame(go_bp@result,
                        stringsAsFactors = F)
    }
    go_cc<-enricher(genelist,
                    pAdjustMethod = 'fdr',minGSSize = 0,
                    # TERM2GENE=PM1GOgene[PM1GOgene$Ontology=='CC',c(2,1)],
                    # TERM2NAME=PM1GOgene[PM1GOgene$Ontology=='CC',c(2,3)],
                    TERM2GENE=PM1GOgene[PM1GOgene$Gene%in%rownames(tclust_cts)&PM1GOgene$Ontology=='CC',c(2,1)],
                    TERM2NAME=PM1GOgene[PM1GOgene$Gene%in%rownames(tclust_cts)&PM1GOgene$Ontology=='CC',c(2,3)],
                    pvalueCutoff=1,
                    qvalueCutoff = 1)
    if(is.null(go_cc)){
      go_cc<-go
    }else{
      go_cc<-data.frame(go_cc@result,
                        stringsAsFactors = F)
    }
    go_mf<-enricher(genelist,
                    pAdjustMethod = 'fdr',minGSSize = 0,
                    # TERM2GENE=PM1GOgene[PM1GOgene$Ontology=='MF',c(2,1)],
                    # TERM2NAME=PM1GOgene[PM1GOgene$Ontology=='MF',c(2,3)],
                    TERM2GENE=PM1GOgene[PM1GOgene$Gene%in%rownames(tclust_cts)&PM1GOgene$Ontology=='MF',c(2,1)],
                    TERM2NAME=PM1GOgene[PM1GOgene$Gene%in%rownames(tclust_cts)&PM1GOgene$Ontology=='MF',c(2,3)],
                    pvalueCutoff=1,
                    qvalueCutoff = 1)
    if(is.null(go_mf)){
      go_mf<-go
    }else{
      go_mf<-data.frame(go_mf@result,
                        stringsAsFactors = F)
    }
    go<-rbind(go_bp,go_cc,go_mf)
    go$cluster<-rep(idx,dim(go)[1])
    return(go)
  }
stopImplicitCluster()
colnames(goresult)
goresult[#goresult$Count>2&
           goresult$p.adjust<0.05&
           goresult$ID%in%slim.go,]%>%
  group_by(cluster)%>%
  # slice_min(order_by = p.adjust, n = 5)%>%
  View(.)