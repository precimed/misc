### This repo contains the scripts to complete the metabolomics projects
### The related documents is:
https://docs.google.com/document/d/1b8a1E55k6djyzuFkl1uU7e8Z6SvJtUMU2aAQ1WtjOjg/edit#heading=h.qjp0rbtewqb
### Running MOSTesT on the UKB data
#### Step 1:
Select phenotypes and create phenotypes csv file by running ukb_helper.py. Example: python ukb_helper.py pheno --out my_pheno --input /ess/p33/data/durable/s3-api/ukblake/phenotypes/*csv \n
--fields 23474 23475 23476 23477 --remove /ess/p33/data/durable/s3-api/ukb_lake/participant_withdrawal/*csv

Create your covariate csv file, e.g. age, sex, 20 genetic principal components: python ukb_helper.py pheno --input /ess/p33/data/durable/s3-api/ukblake/phenotypes/*csv --fields 21022 31 22009 --out my_cov

#### Step 2:
Process your phenotype. Rename columns with actual phenotype names. Make sure to replace any other name that indicates your subjects e.g. "eid" with "IID" and dublicate this column and name it "FID". "FID" and "IID" values has to be integer type. Phenotype csv file must be tab or space delimited. Transform phenotype values by Rank-Based Inverse Normal Transformation.  

#### Step 3:
Create unrelated individuals (king cutoff 0.05) in the genotype file and then extract those induviduals from the imputed files

#### Step 4: 
Run univariate gwas, create permuted genotype (https://github.com/precimed/mostest/blob/mental/mental/permute_bed.py) from the same input (creating chunks of 10K snps using make_chunks_by_snps.py script would be helpful before this)

#### Step 5:
Run univariate gwas on the permuted genotype

#### Step 6:
Extract zscore for all the phenotypes from both univariate and permuted gwas. We will end up with two tables which will be the input for the next step.

#### Step 7:
Run Mostest with the inputs from the previous step.

#### Step 8:
Convert the .mat output into csv using pvals2csv.py script

#### Step 9:
Clump the results with sumstats.py and then annotate the discovery.

#### Nextflow pipeline for steps 3-8:
The nextflow pipeline "mostest_with_nextflow_v1.nf" can be used to complete step 3 to 8 in one go. Read, modify and submit the "implement_nextflow_mostest_v1.job" to accomplish this pipeline. 
