
"""run NMF on sc-RNAseq data.

Usage:
    NMF_scRNAseq.py --min_rank=<int> --max_rank=<int> --max_depth=<int> --RHeatmap=<str> 
        [--method=<str>] [--seed=<str>] [--depth_level=<int>] 
        [--n_run=<int>] [--max_iter=<int>] [--cores=<int>] 
        [--dataSource=<str>] [--algor=<str>]

Options:
    -h --help              # show this screen
    --min_rank <int>       # minimum number of rank
    --max_rank <int>       # maximum number of rank
    --max_depth <int>      # maximum depth to run NMF
    --RHeatmap <str>       # full path to the heatmap R files
    --method <str>         # method implemented in nimfa. lsnmf, bd and nmf (standard NMF) are [default: lsnmf]
    --seed <str>           # seeding choice [default: nndsvd]
    --depth_level <int>    # depth_level, 1 for the start [default: 1]
    --n_run <int>          # number of run [default: 1]
    --max_iter <int>       # maximum number of iterations [default: 50000]
    --cores <int>          # cores to use (not implemented currently) [default: 3]
    --dataSource <str>     # name of the project [default: unknown]
    --algor <str>          # description of running method [default: unknown]

"""

import pandas as pd
import pdb
import subprocess
import glob
import os
import os.path
from docopt import docopt
import nimfa

from contextlib import contextmanager
import timeit
import time

@contextmanager
def time_elapsed(label):
    start_time = timeit.default_timer()
    try:
        yield
    finally:
        elapsed = float("%.2f" % (timeit.default_timer() - start_time))
        # end = time.clock()
        #print "elapsed seconds by %s: %s" %(label, elapsed)
        #print "\n\n"

################################################################################


def FUN_NMF_iterative(
        filename_e,
        filename_o,
        rank,
        nimfa_method=nimfa.Lsnmf, 
        seed="nndsvd",
        depth_level=1, 
        max_depth=1, 
        n_run=1, 
        max_iter=50000 ):

    # pdb.set_trace()
    if max_depth > 0:
        current_depth = max_depth - 1
        # depth_level += 1
        print ("(II) input file:", filename_e)
        ## results will be overwritten!
        ## read the data
        
        npDat= pd.read_csv(filename_e, index_col = 0, sep = ' ', header = 0)
        #npDat=npDat.T
        print(npDat.shape[0])
        ## skip if less than 10 cells exit!
        if npDat.shape[0] < 10:
            print("(II) !Skip as there are less 10 cells in the cluster!")
        else:
            # FUN_estimateRank(npDat, nimfa_method, seed, rank_ranges)
            resWNor=FUN_NMF_iterative_run(npDat, 
                nimfa_method,
                seed, 
                # estimateRank, 
                rank,
                n_run, 
                max_iter,
                filename_o,
                current_depth)
                    # pdb.set_trace()
            FUN_NMF_generate_files(npDat,filename_o,resWNor)
            print ("(II) start new run...")
            FUN_NMF_iterative(
                filename_e,
                filename_o,
                rank,
                nimfa_method,
                seed, 
                depth_level + 1, 
                current_depth, 
                n_run, 
                max_iter)            
        os.chdir("..")
    else:
        print ("(==) reach max depth\n")

################################################################################

def FUN_NMF_iterative_run(npDat, 
        nimfa_method, 
        seed, 
        # estimateRank, 
        rank, 
        n_run, 
        max_iter,
        filename_o,
        current_depth):
    ## process the path to the file;
    # pdb.set_trace()
    dat = npDat
    rowN = list(dat.index)
    colN = dat.columns.values.tolist()

    npDat = dat.values
    n_run = 1 if seed == "nndsvd" else n_run
    # pdb.set_trace()
    #input_prefix = os.path.basename(os.getcwd())
    #fileName = input_prefix + "." + datSource + "." + algor + ".R" + str(rank)
    res = FUN_nmf_run(npDat, nimfa_method, seed, rank, n_run, max_iter)

    NMF_sparW, NMF_sparH = res.fit.sparseness()
    NMF_rss = res.fit.rss()
    NMF_evar = res.fit.evar()
    # print "(==) Algorithm: %s" %(algor)
    #print "(==) rss: %.2f, evar: %.2f, spare_W: %.2f, spare_H: %.2f \n" %(NMF_rss, NMF_evar, NMF_sparW, NMF_sparH)
    ## to start with 1 instead of 0
    metaCells = ["metaCells_" + str(i + 1) for i in range(rank)]

    ## get W matrix
    try:
        resW = pd.DataFrame(res.basis().todense())
        print (" (II) scipy sparse matrixLI_cells")
    except:
        resW = pd.DataFrame(res.basis())
        print ("(II) numpy matrix")
    resW.index = rowN
    resW.columns = metaCells
    #resW.to_csv(fileName + ".basis.csv", sep = "\t", quoting = False) 
    ## normalized by row;
    resWNor = resW.div(resW.sum(axis=1), axis=0)
    #fileName_W = fileName + ".basis.nor.csv"
    #resWNor.to_csv(fileName_W, sep = "\t", quoting = False) 
    

    ## get H matrix
    try:
        resH = pd.DataFrame(res.coef().todense())
        print ("(II) scipy sparse matrix")
    except:
        resH = pd.DataFrame(res.coef())
        print ("(II) numpy matrix")
    resH.columns = colN 
    resH.index = metaCells
    resH = resH.T
    #resH.to_csv(fileName + ".coefficient.csv", sep = "\t", quoting = False) 
    resHNor = resH.div(resH.sum(axis = 1), axis = 0)
    #fileName_H = fileName + ".coefficient.nor.csv"
    #resHNor.to_csv(filename_o, sep = " ", quoting = 0) 
    #resHNor.to_csv('test.csv', sep = " ", quoting = 0)
    if current_depth==0:
        def get_max_index(series):
            #print(series)
            return series.values.argmax()
        resHNor['label']=resHNor.apply(get_max_index,axis=1)
        resHNor.to_csv(filename_o, sep = " ",columns=['label'], quoting = 0) 
    return resWNor

################################################################################
 
"""
subgroups cells with gene expression
"""
def FUN_NMF_generate_files(npDat,filename_o,basis_mat):
    # pdb.set_trace()
    #basisMat_file = glob.glob("*basis.nor.csv")[0]
    #basis_mat = pd.read_csv(basisMat_file, index_col = 0, sep = '\t', header = 0)
    basis_cl = basis_mat.idxmax(axis = 1)
    df_anno = {"cells": basis_cl.index.tolist(), "cls": basis_cl}
    df_anno = pd.DataFrame(df_anno)
    ## convert the basis_cl to a dataFrame in pd
    cls = basis_cl.unique()
    try:
        cls.sort()
    except TypeError:
        pass
    for metaCell in cls:
        LI_cells = df_anno[(df_anno.cls == metaCell)]["cells"].tolist()
        # LI_cells = df_anno[(df_anno.cls == "metaCells_2")]["cells"].tolist()
        ## select cells only in each metacells
        dat_sel = npDat.loc[LI_cells]
        dat_sel = dat_sel[(dat_sel.T != 0).any()]
        ## get the basename of current path
        ## save the file
        #input_prefix = os.path.basename(os.getcwd())
        # pdb.set_trace()
        #dat_sel.to_csv(filename_o, sep = " ", quoting = 0)
    return dat_sel

################################################################################

def FUN_nmf_run(npDat, 
        nimfa_method, 
        seed, 
        rank, 
        n_run, 
        max_iter):
        # **argsToMethod):
    nmfModel = nimfa_method(npDat,
        seed = seed, 
        rank = rank, 
        max_iter = max_iter, 
        n_run = n_run)
        # **argsToMethod)
    #print "(II) rank: %d, method: %s, seed: %s, n_run: %d" %(nmfModel.rank, nmfModel.name, str(nmfModel.seed), nmfModel.n_run)
    print ("(==) starting NMF calculation...")
    nmfFit = nmfModel()
    return nmfFit

################################################################################
#调用格式
#FUN_NMF_iterative('GSE59892_data.csv','result.csv',3)
