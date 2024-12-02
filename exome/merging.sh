#!/bin/bash

# Directory where the files are located
file_directory="/Users/zillurrahman/Desktop/zillur/work/metabolomics_nmr/ukbnmr/burden/results"

# Get a list of unique phenotypes from file names
phenotypes=($(find "$file_directory" -type f -name "ukbnmr_all_gene_burden_chr*_*.regenie" | sed -n 's/.*chr[0-9]*_\(.*\)\.regenie/\1/p' | sort -u))

# Loop through each phenotype and concatenate files, removing the first two lines from all but the first file
for phenotype in "${phenotypes[@]}"; do
  echo $phenotype
  # Concatenate files for the current phenotype and remove the first two lines from the 2nd to the rest
  find "$file_directory" -type f -name "ukbnmr_all_gene_burden_chr*_${phenotype}.regenie" -exec sed '1,2d' {} + | sort -k1,1 -k2,2 | grep -v '^#' > "${phenotype}.csv"
echo 'Done'
done

