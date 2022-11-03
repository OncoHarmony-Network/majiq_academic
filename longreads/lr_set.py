import pandas as pd
from gtfparse import read_gtf

import warnings
warnings.filterwarnings('ignore')

df_gtf = read_gtf('/Users/seongwoohan/Desktop/ex/example.gtf')

junctions = set()
for i in range(len(df_gtf)-1):
    if df_gtf.loc[i,'feature'] == 'exon' and df_gtf.loc[i+1,'feature'] == 'exon':
        #print(df_gtf.loc[i,'feature'], df_gtf.loc[i,'start'], df_gtf.loc[i,'end'])
        junctions.add((df_gtf.loc[i+1,'end'], df_gtf.loc[i,'start']))
        
#print(junctions)
