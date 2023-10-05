### This repo contains the scripts to complete the metabolomics projects
### The related documents is:
https://docs.google.com/document/d/1b8a1E55k6djyzuFkl1uU7e8Z6SvJtUMU2aAQ1WtjOjg/edit#heading=h.qjp0rbtewqb
### Running MOSTesT on the UKB data
#### Step 1:
Select phenotypes and create phenotypes csv file by running ukb_helper.py. Example: python ukb_helper.py pheno --out my_pheno --input /ess/p33/data/durable/s3-api/ukblake/phenotypes/*csv \n
--fields 23474 23475 23476 23477 --remove /ess/p33/data/durable/s3-api/ukb_lake/participant_withdrawal/*csv

Create your covariate csv file, e.g. age, sex, 20 genetic principal components: python ukb_helper.py pheno --input /ess/p33/data/durable/s3-api/ukblake/phenotypes/*csv --fields 21022 31 22009 --out my_cov

#### Step 2:
Process your phenotype. Rename columns with actual phenotype names. Make sure to replace any other name that indicates your subjects e.g. "eid" with "IID" and dublicate this column and name it "FID". "FID" and "IID" values has to be integer type. Phenotype csv file must be tab or space delimited. Transform phenotype values by Rank-Based Inverse Normal Transformation. (If you already have a phenotype and a covariate table ignore these two steps)  

#### Step 3:
Create unrelated individuals (king cutoff 0.05) in the genotype file and then extract those induviduals from the imputed files. In ukb, we merged chr 1-22 together using merge-list in plink and used king 0.05 flag to create a list of unrelated individuals and then extracted those individuals from the imputed bgen files. (If we already have a list of unrelated individuals, this step is not needed) Then we removed rare and duplicated variants (This is mandetory). Example code: 

plink --merge-list ${merge_in} --threads 16 --keep ${pheno} --make-bed --out ${outPrefix}

plink2 --bfile ${inprefix} --king-cutoff 0.05 --make-just-fam --out ${outPrefix} --threads 16 --memory 31000

plink2 --bgen ${bgen} ref-first --sample ${sample} --keep ${fam} --mac 20 --make-bed --out "${chr}_unrelated_imp" --threads 8

plink2 --bed ${beds} --bim ${bims} --fam ${fams} --keep ${pheno} --mac 20 --maf 0.05 --rm-dup 'force-first' --make-bed --out ${chr}_unrelated_nodup --threads 8 --memory 15000

 

#### Step 4: 
Run univariate gwas, create permuted genotype (https://github.com/precimed/mostest/blob/mental/mental/permute_bed.py) from the same input (creating chunks of 10K snps using make_chunks_by_snps.py script would be helpful before this) Example code:
create chunks of 10K snps: 
python make_chunks_by_snps.py ${bims}

plink --bed ${bed} --bim ${bim} --fam ${fam} --extract ${chunk} --make-bed --out "${chunk}" --threads 4

plink2 --bfile ${bed} --glm omit-ref hide-covar --covar ${cov} --covar-variance-standardize --pheno ${pheno} --out ${bed}_glm --threads 4 --memory 7600

python permute_bed.py --bfile ${bed}.csv --out ${bed}_permuted

plink2 --bfile ${bed} --glm omit-ref hide-covar --covar ${cov} --covar-variance-standardize --pheno ${pheno} --out ${bed} --threads 4 --memory 7600


#### Step 5:
Run univariate gwas on the permuted genotype

plink2 --bfile ${bed} --glm omit-ref hide-covar --covar ${cov} --covar-variance-standardize --pheno ${pheno} --out ${bed} --threads 4 --memory 7600

Then merge chunks/chromosome:

python concatenate_chunks.py ${sumstats} ${trait}_combined_original.csv

python concatenate_chunks.py ${sumstats} ${trait}_combined_permuted.csv

##### Upto step 5 is basic. These outputs can used in several downstream analysis 
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
