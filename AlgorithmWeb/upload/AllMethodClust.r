# test SC3
SC3 <- function(filename_e,filename_o,k=1)
{
  library(SingleCellExperiment)
  library(SC3)
  library(scater)
  data <- read.csv(filename_e,sep="")
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
    result <- data.frame(cell_name=colnames(data),cell_label=sce@colData@listData[[2]])
  }
  else
  {
    sce <- sc3(sce, ks = k, biology = FALSE)
    result <- data.frame(cell_name=colnames(data),cell_label=sce@colData@listData[[2]])
  }
  write.csv(result,filename_o,row.names=F)
}



# test CIDR
CIDR <- function(filename_e,filename_o)
{
  library("cidr")
  data <- read.csv(filename_e,sep="")
  sData <- scDataConstructor(as.matrix(data))
  sData <- determineDropoutCandidates(sData)
  sData <- wThreshold(sData)
  sData <- scDissim(sData)
  sData <- scPCA(sData,plotPC=FALSE)
  sData <- nPC(sData)
  sData <- scCluster(sData)
  result <- data.frame(cell_name=colnames(data),cell_label=sData@clusters)
  write.csv(result,filename_o,row.names=F)
}



#test RaceID
RaceID <- function(filename_e,filename_o)
{
  source("/home/projects/AlgorithmWeb/upload/RaceID_class.R")
  x <- read.csv(filename_e,sep="")
  sc <- SCseq(x)
  sc <- filterdata(sc, mintotal=3000, minexpr=1, minnumber=1, maxexpr=500, downsample=FALSE, dsn=1, rseed=17000)
  sc <- clustexp(sc, clustnr=20,bootnr=50,metric="pearson",do.gap=FALSE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=8,rseed=17000)
  sc <- comptsne(sc,rseed=15555)
  sc <- findoutliers(sc, outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.75)
  result <- data.frame(cell_name=names(sc@cpart),cell_label=sc@cpart)
  write.csv(result,filename_o,row.names=F)
}






#test pcaReduce, need k 
#two methods,two results
pcaReduce <- function(filename_e,filename_o1,k,filename_o2="test2.csv")
{
  library("pcaReduce")
  data <- read.csv(filename_e,sep="")
  D <- log2(as.matrix(data) + 1) # log transformation of count data
  Input <- t(D) # data matrix, cells in rows, genes in columns
  Output_S <- PCAreduce(Input, nbt=100, q=30, method='S')
  N <- length(Output_S)
  M <- dim(Output_S[[1]])[2]
  k <- k
  cellL <- Output_S[[N]][,30-k+2]
  result <- data.frame(cell_name=colnames(data),cell_label=cellL)
  write.csv(result,filename_o1,row.names=F)
  
  Output_M <- PCAreduce(Input, nbt=100, q=30, method='M')
  N <- length(Output_M)
  M <- dim(Output_M[[1]])[2]
  k <- k
  cellL <- Output_M[[N]][,30-k+2]
  result <- data.frame(cell_name=colnames(data),cell_label=cellL)
  write.csv(result,filename_o2,row.names=F)
}




# test DIMM_SC, need k, for droplet-base datasets
DIMM_SC <- function(filename_e,filename_o,k)
{
  library("DIMMSC")
  data <- read.csv(filename_e,sep="")
  re <- DIMMSC(data, K=k, method_cluster_initial="kmeans",method_alpha_initial="Ronning", maxiter=200, tol=1e-4, lik.tol=1e-2)
  pre <- re$mem
  result <- data.frame(cell_name=colnames(data),cell_label=pre)
  write.csv(result,filename_o,row.names=F)
}





# test SINCERA,need k
SINCERA <- function(filename_e,filename_o,k)
{
  library("SINCERA")
  expressions <- read.csv(filename_e,sep="")
  data<-expressions
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
  obj <- cluster.geneSelection(obj, method="specificity", pergroup=TRUE, min.samples=min.samples, specifity.thresh=0.7)
  sc <- setGenesForClustering(sc, value=getGenesForClustering(obj))
  sc <- cluster.assignment(sc, k=k, verbose=F)
  result <- data.frame(cell_name=colnames(data),cell_label=sc@data$CLUSTER)
  write.csv(result,filename_o,row.names=F)
}



# test Seurat
Seurat <- function(filename_e,filename_o)
{
  library(Seurat)
  library(dplyr)
  data <- read.csv(filename_e,sep="")
  pbmc <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200)
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
  percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
  pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,do.plot=F,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
  pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 1:5, 
                 genes.print = 5)
  pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
  pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
  pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:20, 
                       resolution = 0.6, print.output = 0, save.SNN = TRUE)
  result <- data.frame(cell_name=colnames(data),cell_label=pbmc@ident)
  write.csv(result,filename_o,row.names=F)
}





# test cellTree,need k
cellTree <- function(filename_e,filename_o,k)
{
  library("cellTree")
  data <- read.csv(filename_e,sep="")
  lda.results <- compute.lda(data,k.topics=k,method = "Gibbs")
  dists <- get.cell.dists(lda.results)
  dists <- as.dist(dists)
  hc <- hclust(dists,method = "complete")
  test <- cutree(hc,k=k)
  result <- data.frame(cell_name=colnames(data),cell_label=test)
  write.csv(result,filename_o,row.names=F)
}





#test SNN_Clique
SNN <- function(filename_e,filename_o="edge.txt")
{
  source('/home/projects/AlgorithmWeb/upload/SNN.R')
  data <- read.csv(filename_e,sep="")
  data<-log2(data+1)
  data <- t(data)
  SNN(data,filename_o, k=3, distance='euclidean')
}
# perform snncliq.py
# cd /Users/ruiyi/Documents/Project/CLUST/code/TestNormalization
# python Cliq.py -i edge.txt -o output.txt
output <- function(filename_e,filename_i="output.txt",filename_o)
{
  data <- read.csv(filename_e,sep="")
  p_clust <- read.csv(filename_i, sep="",header=F)
  max_value <- max(p_clust)
  ln <- length(which(p_clust[,1]==-1))
  if(ln>0)
  {
    p_clust[which(p_clust[,1]==-1),] <- c((max_value+1):(max_value+ln))
  }
  result <- data.frame(cell_name=colnames(data),cell_label=p_clust$V1)
  write.csv(result,filename_o,row.names=F)
  #adjustedRandIndex(o_clust$clustLable,p_clust$V1)
}




#test giniclust
giniClust <- function(filename_e,filename_o)
{
  data.file = filename_e
  exprimentID = basename(data.file)
  source("/home/projects/AlgorithmWeb/upload/GiniClust_parameters.R")
  source("/home/projects/AlgorithmWeb/upload/GiniClust_packages.R")
  #preprocess
  source("/home/projects/AlgorithmWeb/upload/GiniClust_Preprocess.R")
  ExprM.Results = GiniClust_Preprocess(data.file,exprimentID)
  ExprM.RawCounts = ExprM.Results$raw
  #filtering
  source("/home/projects/AlgorithmWeb/upload/GiniClust_Filtering.R")
  data.file = filename_e
  exprimentID = basename(data.file)
  ExprM.Results.filter = GiniClust_Filtering(ExprM.RawCounts,exprimentID)
  ExprM.RawCounts.filter = ExprM.Results.filter$raw
  #gene selection
  source("/home/projects/AlgorithmWeb/upload/GiniClust_Fitting.R")
  GeneList.final = GiniClust_Fitting(ExprM.RawCounts.filter,exprimentID)
  #clustering
  source("/home/projects/AlgorithmWeb/upload/GiniClust_Clustering.R")
  Cluster.Results = GiniClust_Clustering(ExprM.RawCounts,ExprM.RawCounts.filter,GeneList.final,eps,MinPts,exprimentID)
  cell.cell.distance = Cluster.Results$cell_cell_dist
  c_membership = Cluster.Results$c_membership
  clustering_membership_r = Cluster.Results$clustering_membership_r
  rare.cells.list.all = Cluster.Results$rare.cell
  filename = paste(exprimentID,"_clusterID.csv",sep='')
  result <- read.csv(filename)
  result <- data.frame(cell_name=result$cell.ID,cell_label=result$cluster.ID)
  write.csv(result,filename_o,row.names=F)
}

BISCUIT <- function(filename_e,filename_o)
{
  library(MCMCpack)
  library(mvtnorm)
  library(ellipse)
  library(coda)
  library(Matrix)
  library(Rtsne)
  library(gtools)
  library(foreach)
  library(doParallel)
  library(doSNOW)
  library(snow)
  library(lattice)
  library(MASS)
  library(bayesm)
  library(robustbase)
  library(chron)
  library(mnormt)
  library(schoolmath)
  library(RColorBrewer)
  cat('file-name-is')
  cat(filename_e)
  input_file_name<- filename_e
  data <- read.csv(filename_e,sep="")
  input_data_tab_delimited <- TRUE #set to TRUE if the input data is tab-delimited
  is_format_genes_cells <-  TRUE #set to TRUE if input data has rows as genes and columns as cells
  choose_cells <- ncol(data) #comment if you want all the cells to be considered
  choose_genes <- 150 #comment if you want all the genes to be considered
  gene_batch <- 50 #number of genes per batch, therefore num_batches = choose_genes (or numgenes)/gene_batch. Max value is 150
  num_iter <- 20 #number of iterations, choose based on data size.
  num_cores <- detectCores() - 4 #number of cores for parallel processing. Ensure that detectCores() > 1 for parallel processing to work, else set num_cores to 1.
  z_true_labels_avl <- FALSE #set this to TRUE if the true labels of cells are available, else set it to FALSE. If TRUE, ensure to populate 'z_true' with the true labels in 'BISCUIT_process_data.R'
  num_cells_batch <- ncol(data) #set this to 1000 if input number of cells is in the 1000s, else set it to 100.
  alpha <- 0.005 #DPMM dispersion parameter. A higher value spins more clusters whereas a lower value spins lesser clusters.
  output_folder_name <- "output" #give a name for your output folder.
  cat(filename_e)
  source("/home/projects/AlgorithmWeb/upload/BISCUIT/BISCUIT_main.R")
  cat(filename_e) 
  #compute result
  filename <- paste(output_folder_name,"/home/projects/AlgorithmWeb/cluster_probabilities.csv",sep='')
  p_clust <- read.csv(filename,row.names=1)
  result <- data.frame(cell_name=colnames(data),cell_label=p_clust$z_inferred)
  write.csv(result,filename_o,row.names=F)
}
