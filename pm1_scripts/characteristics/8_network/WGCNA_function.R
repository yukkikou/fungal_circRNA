## function
# used in PCA_cluster.R and WGCNA.R
callog2Cpm = function(mtx){
  mtx_tmp = mtx[, -1]
  libSize = colSums(mtx_tmp)
  mtx_out = data.frame(matrix(nrow = dim(mtx)[1], ncol = dim(mtx)[2]))
  for (i in 1:dim(mtx_tmp)[2]){
    mtx_out[, (i+1)] = log2(mtx_tmp[, i]/libSize[i]*1000000+1)
  }
  print(colSums(mtx_out))
  mtx_out[, 1] = mtx[, 1]
  colnames(mtx_out) = colnames(mtx)
  as_tibble(mtx_out)
}

calCpm = function(mtx){
  mtx_tmp = mtx[, -1]
  libSize = colSums(mtx_tmp)
  mtx_out = data.frame(matrix(nrow = dim(mtx)[1], ncol = dim(mtx)[2]))
  for (i in 1:dim(mtx_tmp)[2]){
    mtx_out[, (i+1)] = mtx_tmp[, i]/libSize[i]*1000000
  }
  print(colSums(mtx_out))
  mtx_out[, 1] = mtx[, 1]
  colnames(mtx_out) = colnames(mtx)
  as_tibble(mtx_out)
}

mtxFilter = function(mtx, thd, colData, prec){
  mtx_tmp = mtx %>%
    pivot_longer(cols = colnames(mtx)[-1], 
                 names_to = "sample",
                 values_to = "count") %>%
    left_join(colData, by = c("sample" = "sample")) %>%
    mutate(keep = as.numeric(count > thd)) %>%
    group_by(feat_id, phase) %>%
    summarise(k_score = sum(keep)/3,
              k_keep = as.numeric(k_score >= prec)) %>%
    group_by(feat_id) %>%
    summarise(f_score = sum(k_keep),
              f_keep = (f_score) > 0) %>%
    filter(f_keep == T) %>%
    as.data.frame()
  
  mtx %>% filter(feat_id %in% mtx_tmp$feat_id)
}

top1000zs = function(mtx, thd, colData, prec, number){
  mtx_tmp = mtx
  thd_tmp = thd
  colData_tmp = colData
  prec_tmp = prec
  # expression normalization  
  cpm = calCpm(mtx)
  
  # expression filtering and low variance
  cpmFvars = mtxFilter(mtx = mtx_tmp, thd = thd_tmp, colData = colData_tmp, prec = prec_tmp) %>%
    column_to_rownames("feat_id") %>%
    apply(., 1, var)
  
  keep_feat = names(sort(cpmFvars, decreasing = T)[1:number])
  
  # z-score
  keep_cpm = cpm %>% as.data.frame() %>%
    filter(feat_id %in% keep_feat) 
  
  keep_cpm %>%
    column_to_rownames("feat_id") %>%
    t() %>%
    apply(., 1, function(x) scale(x , center = T, scale = T)) %>%
    `rownames<-`(keep_cpm$feat_id)
}

pcaPlot = function(id_mtx, prefix, colData, leading_number = 100){
  PC = prcomp(t(id_mtx))
  message(paste0("Generate variable: ", prefix, "_pca"))
  assign(paste0(prefix, "_pca"), PC, envir = .GlobalEnv)
  
  PC_sd <- setNames(PC$sdev , paste0("PC",1:length(PC$sdev)))
  PC_var_expl <- round(( PC_sd^2 ) / sum(PC_sd^2) * 100, 2)
  
  xlab = paste0("PC1 ","(", paste0(as.character(PC_var_expl[1]), "%", ")", sep=""))
  ylab = paste0("PC2 ","(", paste0(as.character(PC_var_expl[2]), "%", ")", sep=""))
  if (grepl(prefix, pattern = "circ")){
    cols = col_nejm[1:3]
  } else {
    cols = col_nejm[4:6]
  }
  
  plot_pca = PC$x[,1:2] %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(colData) %>%
    ggplot(aes(x = PC1, y = PC2, color = label, label = sample)) +
    geom_point(size = 1.5) +
    geom_text(nudge_x = 1, nudge_y = 1, check_overlap = T, size = 2.5) +
    labs(x = xlab, y = ylab) +
    scale_color_manual(values = cols) +
    theme_bw() +
    mytheme +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )
  
  message(paste0("Generate variable: ", "plot_", prefix, "_pca"))
  assign(paste0("plot_", prefix, "_pca"), plot_pca, envir = .GlobalEnv)
  
  leading_genes = PC$rotation
  
  leading_PC1 = sort(leading_genes[,1],decreasing = T)[1:leading_number] %>%
    as.data.frame() %>%
    rownames_to_column("feat_id") %>%
    `colnames<-`(c("feat_id", "rotation")) %>%
    mutate(pc = "PC1")
  leading_PC2 = sort(leading_genes[,2],decreasing = T)[1:leading_number] %>%
    as.data.frame() %>%
    rownames_to_column("feat_id") %>%
    `colnames<-`(c("feat_id", "rotation")) %>%
    mutate(pc = "PC2")
  
  leading_PC = rbind(leading_PC1, leading_PC2) %>% as_tibble()
  message(paste0("Generate variable: ", prefix, "_pca_leadingGene"))
  assign(paste0(prefix, "_pca_leadingGene"), leading_PC, envir = .GlobalEnv)
}


datExprClean = function(expr = datExpr0){
  #check missing and outlier
  gsg = goodSamplesGenes(expr, verbose = 3)
  gsg$allOK
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(expr)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(expr)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    expr = expr[gsg$goodSamples, gsg$goodGenes]
  }
  #cluster the samples to see if there are any obvious outliers
  sampleTree = hclust(dist(expr), method = "average")
  # sizeGrWindow(12,9);par(cex = 0.6);par(mar = c(0,4,2,0))
  # plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
  #      cex.axis = 1.5, cex.main = 2)
  # abline(h = 100, col = "red")
  # clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
  # table(clust)
  # keepSamples = (clust==1)
  # datExpr = datExpr0[keepSamples, ]
  nGenes <<- ncol(expr)
  nSamples <<- nrow(expr)
  collectGarbage()
  expr
}

sampleTree = function(expr = datExpr, traitColors, anno){
  sampleTree2 = hclust(dist(expr), method = "average")
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = anno)
}

powerSelect = function(expr = datExpr){
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5)
  # Plot the results:
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2))
  cex1 = 0.9
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  message(paste("Suggested softthd is", sft$powerEstimate))
  
  sft$powerEstimate
}

dyMotif = function(expr = datExpr, softPower = sft){
  adjacency = adjacency(expr, power = softPower)
  
  TOM <<- TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  
  geneTree <<- hclust(as.dist(dissTOM), method = "average")
  
  sizeGrWindow(12,9)
  plot(geneTree, xlab="", sub="", 
       main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)
  
  minModuleSize = 30
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  table(dynamicMods)
  
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  
  # Calculate eigengenes
  MEList = moduleEigengenes(expr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  # Plot the result
  sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  
  #as figure shows
  MEDissThres = 0.15
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(expr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  sizeGrWindow(12, 9)
  #pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  #dev.off()
  # Rename to moduleColors
  
  moduleColors <<- mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50))
  moduleLabels <<- match(moduleColors, colorOrder)-1
  MEs <<- mergedMEs
  tmp = list(moduleColors, moduleLabels, MEs, TOM)
  names(tmp) = c("moduleColors", "moduleLabels", "MEs", "TOM")
  tmp
}

toColData = function(expr = datExpr, mcol = moduleColors, ns = nSamples, anno = colDataMtx){
  MEs0 = moduleEigengenes(expr, mcol)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, anno, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, ns)
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(anno),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  moduleTraitCor_tmp = moduleTraitCor %>%
    as.data.frame() %>%
    rownames_to_column("module")
  
  moduleTraitPvalue_tmp = moduleTraitPvalue %>%
    as.data.frame() %>%
    rownames_to_column("module")
  
  left_join(moduleTraitCor_tmp, moduleTraitPvalue_tmp,
            by = c("module" = "module"),
            suffix = c(".cor", ".pvalue"))
}

calGeneModuleColdata = function(mes = MEs, expr = datExpr, ns = nSamples, anno, prefix, moduleCol){
  modNames = substring(names(mes), 3)
  geneModuleMembership = as.data.frame(cor(expr, mes, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), ns))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  
  # Define continuous variable
  pheno = as.data.frame(anno)
  names(pheno) = prefix
  
  geneTraitSignificance = as.data.frame(cor(expr, pheno, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), ns))
  names(geneTraitSignificance) = paste("GS.", prefix, sep="")
  names(GSPvalue) = paste("p.GS.", prefix, sep="")
  
  #take module brown as example, identifying genes with high GS and MM
  module = moduleCol
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  sizeGrWindow(7, 7)
  par(mfrow = c(1,1))
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for body weight",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)  
}

exGeneFromModule = function(expr = datExpr, mod_cols = moduleColors, moduleCol){
  probes = colnames(expr)
  inModule = (mod_cols==moduleCol)
  probes[inModule]
}

exGeneModuleCol = function(mod_cols = moduleColors, expr = datExpr, search = T, feat){
  if (search){
    print(grep(colnames(expr), pattern = feat, value = T))
    mod_cols[grep(colnames(expr), pattern = feat)]
  } else {
    mod_cols[(colnames(expr)== feat)]
  }
}

outModuleToCyto = function(expr = datExpr, tom = TOM, mod_cols = moduleColors, moduleCol, filter = 30){
  probes = colnames(expr)
  inModule = (mod_cols==moduleCol)
  modProbes = probes[inModule]
  
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  
  if (filter){
    nTop = filter
    IMConn = softConnectivity(expr[, modProbes])
    top = (rank(-IMConn) <= nTop)
    modTOM = modTOM[top, top]
  }
  
  cyt = exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(moduleCol, collapse="-"), ".txt", sep=""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(moduleCol, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = modProbes, 
    nodeAttr = moduleColors[inModule]
  )
}
