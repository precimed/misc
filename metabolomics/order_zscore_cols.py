import pandas as pd
import sys

# Get the file paths from command line arguments
file1 = sys.argv[1]
file2 = sys.argv[2]
output1 = sys.argv[3]
output2 = sys.argv[4]

# Read the CSV files into DataFrames
df1 = pd.read_csv(file1, sep='\t', index_col = False)
df2 = pd.read_csv(file2, sep='\t', index_col = False)

# Get the ordered column names from df1
ordered_columns = df1.columns.tolist()

# Reorder columns in df2_renamed to match the order in df1_renamed
df2_reordered = df2[ordered_columns]

# Write the ordered data frames to CSV files
df1=df1.drop(0)
df1.to_csv(output1, index=False, sep='\t')
df2_reordered=df2_reordered.drop(0)
df2_reordered.to_csv(output2, index=False, sep='\t')
