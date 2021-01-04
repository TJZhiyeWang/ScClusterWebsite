RG <- function(filename_e,k,filename_o)
{
  data <- read.csv(filename_e,sep="")
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if(!requireNamespace("SC3", quietly = TRUE))
  {
    BiocManager::install("SC3")
  }
  if(!requireNamespace("SingleCellExperiment", quietly = TRUE))
  {
    BiocManager::install("SingleCellExperiment")
  }
  if(!requireNamespace("scater", quietly = TRUE))
  {
    BiocManager::install("scater")
  }
  library(SingleCellExperiment)
  library(SC3)
  library(scater)
  #fake cell types
  ann <- data.frame(cell_type1=c(1:ncol(data)))
  row.names(ann) <- colnames(data)
  sce <- SingleCellExperiment(
    assays = list(
      counts = as.matrix(data),
      logcounts = log2(as.matrix(data) + 1)
    ),
    colData = ann
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)
  if(k == 1)
  {
    sce <- sc3(sce, k_estimator = TRUE, biology = FALSE)
    r_sc3 <- data.frame(cell_name=colnames(data),cell_label=sce@colData@listData[[2]])
  } else
  {
    sce <- sc3(sce, ks = k, biology = FALSE)
    r_sc3 <- data.frame(cell_name=colnames(data),cell_label=sce@colData@listData[[2]])
  }
  
  
  
  if (!suppressWarnings(require("cidr",quietly = TRUE))) {
    install.packages("cidr",dependencies = TRUE, repos="http://cran.r-project.org")
  }
  library(cidr,quietly = TRUE)
  sData <- scDataConstructor(as.matrix(data))
  sData <- determineDropoutCandidates(sData)
  sData <- wThreshold(sData)
  sData <- scDissim(sData)
  sData <- scPCA(sData,plotPC=FALSE)
  sData <- nPC(sData)
  sData <- scCluster(sData)
  r_cidr <- data.frame(cell_name=colnames(data),cell_label=sData@clusters)
  
  
  if (!suppressWarnings(require("SINCERA",quietly = TRUE))) {
    install.packages("SINCERA",dependencies = TRUE, repos="http://cran.r-project.org")
  }
  library(SINCERA,quietly = TRUE)
  expressions <- data
  cells <- list(CELL=colnames(expressions),SAMPLE=matrix(1,ncol(expressions),1),CLUSTER=c(1:ncol(data)))
  genes <- list(SYMBOL=rownames(expressions))
  sc <- construct(exprmatrix=expressions,
                  samplevector=paste("sample", cells$SAMPLE, sep=""))
  sc <- prefilterGenes(sc, pergroup=TRUE, min.expression=5, min.cells=2, min.samples=1)
  sc <- expr.minimum(sc, value=0.01)
  sc <- normalization.zscore(sc, pergroup=TRUE)
  sc <- doPCA(sc, genes=NULL, use.fast = T)
  sc <- doTSNE(sc, genes=NULL, dims = 1:5,perplexity=10, use.fast = T)
  min.samples <- 1
  obj <- prefilterGenes(sc, pergroup=TRUE, min.expression=1, min.cells=10, min.samples=min.samples)
  obj <- cluster.geneSelection(obj, method="specificity", pergroup=TRUE, min.samples=min.samples, specifity.thresh=0.7,do.plot = F)
  # genes with at least 0.7 specificity in at least 6 samples
  obj <- cluster.geneSelection(obj, method="specificity", pergroup=TRUE, min.samples=min.samples, specifity.thresh=0.7,do.plot=F)
  sc <- setGenesForClustering(sc, value=getGenesForClustering(obj))
  sc <- cluster.assignment(sc, k=k, verbose=F)
  r_sincera <- data.frame(cell_name=colnames(data),cell_label=sc@data$CLUSTER)
  
  
  library("pcaReduce")
  D <- log2(as.matrix(data) + 1) # log transformation of count data
  Input <- t(D) # data matrix, cells in rows, genes in columns
  Output_S <- PCAreduce(Input, nbt=1, q=45, method='S')
  N <- length(Output_S)
  M <- dim(Output_S[[1]])[2]
  k <- k
  cellL <- Output_S[[N]][,45-k+2]
  r_pcaReduce <- data.frame(cell_name=colnames(data),cell_label=cellL)
  

  cn <- nrow(r_sc3)
  sim_M <- matrix(0,cn,cn)
  row.names(sim_M) <- r_sc3$cell_name
  colnames(sim_M) <- r_sc3$cell_name
  allResult <- list(r1=r_sc3,r2=r_cidr,r3=r_sincera,r4=r_pcaReduce)
  Cells <- as.matrix(row.names(sim_M))
  for (m in 1:4)
  {
    temp <- allResult[[m]]
    ln <- length(unique(temp$cell_label))
    for (c in 1:ln)
    {
      cell <- as.matrix(temp$cell_name[which(temp$cell_label==c)])
      loc <- matrix(0,nrow(cell),1)
      for (i in 1:nrow(cell))
      {
        loc[i,1] <- which(Cells==cell[i,1])
      }
      sim_M[loc,loc] <- sim_M[loc,loc]+1
    }
  }
  test <- sim_M
  test <- as.matrix(test)
  for (i in 1:cn)
  {
    test[which(test[,i]<2),i] <- 0
  }
  
  if (!suppressWarnings(require("igraph",quietly = TRUE))) {
    install.packages("igraph",dependencies = TRUE, repos="http://cran.r-project.org")
  }
  library(igraph,quietly = TRUE)
  g <- graph_from_adjacency_matrix(test, weighted=TRUE,mode="undirected")
  re2<-walktrap.community(g,weights=E(g)$weight,step=4)
  result <- data.frame(cell_name=colnames(data),cell_label=re2$membership)
  write.csv(result,filename_o,row.names=F)
  re <- list(sim_M=sim_M,r1=r_sincera,r2=r_cidr,r3=r_sc3,r4=r_pcaReduce,r5=result)
  return(re)
}
