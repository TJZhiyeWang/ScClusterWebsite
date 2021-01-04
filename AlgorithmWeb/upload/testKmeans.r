
KMeans<-function(filename_e,c_number,filename_output){

    data <- read.table(filename_e,sep="")
    data <- as.matrix(data)
    data <- t(data)
    re_clust <- kmeans(x=data, centers=c_number, iter.max = 1000, nstart = 10,algorithm = "Hartigan-Wong", trace=FALSE)
    cell_label <- re_clust$cluster
    result <- data.frame(cell_name=rownames(data),cell_label=cell_label)
    write.csv(result,filename_output)
}
