# Authors: Lan Jiang; Huidong Chen
# Contact information: lan_jiang@hms.harvard.edu; 
#                      hdchen@jimmy.harvard.edu;

#########################################################
### All the parameters used in Giniclust
#########################################################
data.type = 'RNA-seq'
if(data.type == 'RNA-seq'){
  minCellNum           = 3                                                # filtering, for at least expressed in how many cells
  minGeneNum           = 2000                                             # filtering, for at least expressed in how many genes
  expressed_cutoff     = 0                                                # filtering, for raw counts
  log2.expr.cutoffl    = 0                                                # cutoff for range of gene expression   
  log2.expr.cutoffh    = 20                                               # cutoff for range of gene expression 
  Gini.pvalue_cutoff   = 0.01                                           # fiting, Pvalue, control how many gene finally used.
  Norm.Gini.cutoff     = 1                                                # fiting, NomGini, control how many gene finally used, 1 means not used.
  span                 = 0.9                                              # parameter for LOESS fitting
  outlier_remove       = 0.75                                             # parameter for LOESS fitting
  Gamma                = 0.9                                              # parameter for clustering
  diff.cutoff          = 1                                                # MAST analysis, filter gene don't have high log2_foldchange to reduce gene num
  lr.p_value_cutoff    = 1e-5                                             # MAST analysis, pvalue cutoff to identify differential expressed gene
  CountsForNormalized  = 100000                                           
  rare_p               = 0.05                                             # propostion of cell number < this value will be considered as rare cell clusters.
  perplexity           = 30
  eps                  = 0.5
  MinPts               = 3
}
