import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial

chunk_size=1000000
#for i in range(np.round(2167207740/chunk_size).astype('int')):

def take_uniq(i,inputA,outputA):
    data=pd.read_csv(inputA,sep=',',names=['only'],header=None,nrows=chunk_size,skiprows=chunk_size*i)
    data = data.sort_values(by='only')
    first = data.groupby('only').first().reset_index()
    print(i/np.round(2167207740/chunk_size))
    first.to_csv(outputA,sep='\t',index=False,header=False,mode='a')

indices=range(np.round(2167207740/chunk_size).astype('int'))

pool = mp.Pool(20)
res  = pool.map(partial(take_uniq,inputA='slop_mega_BS_Ch1.txt',outputA='slop_uniq_mega_BS_Ch1.txt'),indices)

