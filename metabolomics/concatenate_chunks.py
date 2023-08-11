import pandas as pd
import sys
import os

# Read the file paths and output file from command line arguments
file_paths = sys.argv[1:-1]
output_file = sys.argv[-1]

dfs=[]
for f in file_paths:
    df=pd.read_csv(f, index_col=False, sep='\t') 
    if len(df) >= 2:
        dfs.append(df)

combined_df = pd.concat(dfs, ignore_index=True)

# Sort the DataFrame by the first and second columns
sorted_df = combined_df.sort_values(by=[combined_df.columns[0],combined_df.columns[1]])

# Write the sorted DataFrame to the output file
sorted_df.to_csv(output_file, index=False,sep='\t')
print(output_file, " Done!")
