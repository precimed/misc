# misc
Misc of analysis

### Running MOSTesT for UKB data
#### Step 1:
Select phenotype, run ukb_helper.py. Example: python ukb_helper.py pheno --out my_pheno --input /cluster/projects/p33/s3-api/ukblake/phenotypes/*csv \n
--fields 31 22001 22006 21022 30620 30600 --remove /cluster/projects/p33/s3-api/ukb_lake/participant_withdrawal/*csv


#### Step 2:
Process your phenotype. Rename columns with actual phenotype names. Transform phenotype values.

#### Step 3:
Create unrelated individuals in the genotype file and then extract those induviduals from the imputed files

#### Step 4: 
Run univariate gwas, create permuted genotype from the same input

#### Step 5:
Run univariate gwas on the permuted genotype

#### Step 6:
Dark....
