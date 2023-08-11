import pandas as pd
import os
import sys

# Read file names from command line arguments
file_names = sys.argv[1:-1]
output_file_path = sys.argv[-1]

print("The files are:", len(file_names))

# Initialize an empty list to store the data frames
data_frames = []

# Process each file
for file_name in file_names:
    # Read the 11th column of the file
    data = pd.read_csv(file_name, usecols=[10], header=None, sep="\t")
    print("Processing:", file_name)
    data_frames.append(data)

# Concatenate the data frames
combined_data = pd.concat(data_frames, axis=1)

print("The size of data:", combined_data.shape)
print(combined_data.head())

# Set column names
col_names = [f.replace('_combined_original.csv','').replace('_combined_permuted.csv', '') for f in file_names]
#with open(output_file_path, 'w') as outf:
#    outf.write('\t'.join(col_names) + '\n')
print("Columns of the data:", col_names)
combined_data.columns = col_names

# Write output file
combined_data.to_csv(output_file_path, sep="\t", index=False)
