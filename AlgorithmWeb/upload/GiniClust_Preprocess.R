# Author: Lan Jiang
# Contact information: lan_jiang@hms.harvard.edu

#################### Preprocessing ####################
GiniClust_Preprocess <- function(data.file,exprimentID){
  data.type = 'RNA-seq'
  if(data.type == 'RNA-seq'){
    ExprM.RawCounts  <- read.csv(data.file, sep=" ", head=T)
    #write.table(ExprM.RawCounts, file=paste(out.folder, "/", exprimentID, "_rawCounts.csv", sep=""), sep=",", row.names = TRUE, col.names = TRUE, quote=FALSE)
    return(list('raw' = ExprM.RawCounts))
  }
}
####################      End     ####################
