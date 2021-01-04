import numpy as np
import pandas as pd
from backspinpy import backSPIN
def BACKSPIN(filename_e,filename_o,kvalue):
    df = pd.read_csv(filename_e,sep=' ')
    data=df.values
    ans = backSPIN(data, numLevels=kvalue, verbose=True)
    seq = ans.cells_bor_level[len(ans.cells_bor_level)-1]
    ind = ans.cells_order
    cellname=df.columns.to_list()
    label=[]
    cell=[]
    label_ind=0
    for i in range(seq.shape[0]-1):
        n=seq[i]
        while n<seq[i+1]:
            label.append(label_ind)
            cell.append(cellname[n])
            n=n+1
        label_ind+=1
    result=pd.DataFrame({'cell':cell,'label':label})
    result.set_index('cell',inplace=True)
    result.to_csv(filename_o,sep=' ')
