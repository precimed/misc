import pandas as pd
import numpy as np
import sys
import os
import glob


data=sys.argv[1]
df=pd.read_csv(data,sep='\t',header=None)
df=df.iloc[:, 1]
st1=data.split('_')[0]

max_rows = 10000
dataframes = []
while len(df) > max_rows:
    top = df[:max_rows]
    dataframes.append(top)
    df = df[max_rows:]
    #print(df.head())
else:
    dataframes.append(df)
#You could then save out these data frames:

for _, frame in enumerate(dataframes):
    frame.to_csv(f'{st1}'+'_'+'chunk'+str(_)+'.csv', index=False,sep='\t',header=False)
