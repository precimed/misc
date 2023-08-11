import pandas as pd
import numpy as np
import sys
import os

data=sys.argv[1]
df=pd.read_csv(data,sep='\t',header=None)
chnk_ult=df.iloc[:,3].max()
bim_arr=df.iloc[:,3].values
chnk_size=300000
chnk_start=np.arange(0,chnk_ult,chnk_size)
chnk_end=chnk_start+chnk_size

df2=pd.DataFrame()
df2['x']=chnk_start+1
df2['y']=chnk_end
df2[f'{data}']=os.path.basename(data)
df2.iloc[:,2]=df2.iloc[:,2].str.split('.').str[0]
df2=df2[((df2['x'].values[:, None] <= bim_arr) & (df2['y'].values[:, None] >= bim_arr)).any(1)]
df2.index.name = 'chunk_array'
df2.reset_index(inplace=True)
order = [3,0,1,2]
df2 = df2[[df2.columns[i] for i in order]]
df2 = df2.replace({'chr':'','_unrelated_imp':''})
df2.iloc[:,0]=df2.iloc[:,0].str.replace('chr','')
df2.iloc[:,0]=df2.iloc[:,0].str.replace('_unrelated_imp','')
df3=df2.iloc[:,[0,1,2,3]]
#print(df2)
df3.to_csv(f'{data}_chunks_for_permute.csv',index=False,header=False,sep='\t')
