## 24th May 2017
## BISCUIT postprocessing
##
## Code author SP
##

###
###

if(num_gene_batches ==1){
    
    source("BISCUIT_Rfunction/BISCUIT_parallel_impute_onegenebatch.R")
    source("BISCUIT_Rfunction/BISCUIT_extras_onegenebatch.R")
 
}else{
    
    source("BISCUIT_Rfunction/BISCUIT_parallel_impute.R")
    source("BISCUIT_Rfunction/BISCUIT_extras.R")

}

