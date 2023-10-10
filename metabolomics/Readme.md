### This repo contains an overview of the steps (incl. required code) to run uni- and multivariate GWAS. 
Here, we focus on Nightingale NMR data, though of course the steps involved are the same for any set of measures suitable for multivariate GWAS

#### Step 1: Select phenotype and covariate data
Select all 249 Nightingale markers and create a 'phenotypes' dataframe (rows subjects, columns markers).
Additionally, create a covariates dataframe. In this project, we use age, sex, and the first 20 genetic principal components.

For UKB data, we would use the [ukb_helper.py](https://github.com/precimed/ukb) tool, as follows:

```
python ukb_helper.py pheno --out my_pheno --input /ess/p33/data/durable/s3-api/ukblake/phenotypes/*csv \n
--fields 23474 23475 23476 23477 --remove /ess/p33/data/durable/s3-api/ukb_lake/participant_withdrawal/*csv

python ukb_helper.py pheno --input /ess/p33/data/durable/s3-api/ukblake/phenotypes/*csv --fields 21022 31 22009 --out my_cov
```

Note: in this project, we are interested in two variations on the main analysis: 1) sex-specific, and 2) covaried for BMI. 

For 1): simply create two phenotype files, one for each sex, and run through the steps below for both of them. Drop the sex covariate from covariates file

For 2): run through the exact same steps as for the main analysis, but include BMI in the covariates file. In UKB field 21001.

#### Step 2: Some light formatting of this data
Ensure phenotype columns have informative names. For UKB, we replace field IDs with biomarker names at this point (underscores okay, no other special characters). See the 'dictionary' file for the names we use. In line with Plink, subject ID columns in both the phenotype and covariate dataframe should be named "IID". Duplicate this column and name it "FID".
Apply inverse rank-based normal transformation to all the phenotype columns. In R, the code for this (whereby df is your phenotype dataframe, and features is the vector of marker/column names) is:

```
for(f in features){df[,f] <- qnorm(ecdf(df[,f])(df[,f]) - 0.5/dim(df)[1])}
```

Write the resulting data to text files (one for phenotypes, one for covariates), tab or space delimited.

#### Step 3: Readying genotype data
Create a list of unrelated individuals (king cutoff 0.05) based on the genotype file and then extract those individuals from the imputed files. In UKB, we first merged chr 1-22 together.  Example code: 

```
plink --merge-list ${merge_in} --threads 16 --keep ${pheno} --make-bed --out ${merged_geno}
plink2 --bfile ${merged_geno} --king-cutoff 0.05 --make-just-fam --out ${unrelated_geno} --threads 16 --memory 31000
plink2 --bgen ${bgen} ref-first --sample ${sample} --keep ${unrelated_geno} --make-bed --out ${chr}_unrelated_imp --threads 8
```

Apply some variant filtering, removing variants that are rare (maf .005), duplicated, or often missing (geno .05). We also propose to include a filter for removing individuals with much missingness (mind .1). We do note that applicability of such filters is a bit dependent on sample characteristics, data quality, and previous QC steps. 

```
plink2 --bfile ${chr}_unrelated_imp --keep ${pheno} --maf 0.005 --rm-dup 'force-first' --geno .05 --mind .1 --make-bed --out ${chr}_unrelated_nodup --threads 8 --memory 15000
```

#### Step 4: Run univariate GWAS on permuted and unpermuted genotype data
To aid computation/permutation, first create chunks of 10K snps using [make_chunks_by_snps.py](https://github.com/precimed/misc/blob/main/metabolomics/make_chunks_by_snps.py) script:  

```
python make_chunks_by_snps.py ${chr}_unrelated_imp_nodup.bim
```

Run univariate GWAS on (chunks of) unpermuted data:

```
plink --bfiles ${chr}_unrelated_imp_nodup --extract ${chunk} --make-bed --out "${chunk}" --threads 4
plink2 --bfile ${chunk} --glm omit-ref hide-covar --covar ${cov} --covar-variance-standardize --pheno ${pheno} --out ${bed}_glm --threads 4 --memory 7600
```

Create permuted genotype through our [permute_bed](https://github.com/precimed/mostest/blob/mental/mental/permute_bed.py) tool and run permuted GWAS:

```
python permute_bed.py --bfile ${chunk} --out ${chunk}_permuted
plink2 --bfile ${chunk}_permuted --glm omit-ref hide-covar --covar ${cov} --covar-variance-standardize --pheno ${pheno} --out ${bed} --threads 4 --memory 7600
```

#### Step 5: Merge chunks/chromosome:
apply our [concatenate tool](https://github.com/precimed/misc/blob/main/metabolomics/concatenate_chunks.py): 

```
python concatenate_chunks.py ${sumstats} ${trait}_combined_original.csv
python concatenate_chunks.py ${sumstats} ${trait}_combined_permuted.csv
```

##### That was it! At least for the purposes of creating the sumstats files needed for downstream analyses. Below steps for running MOSTest (to do: write up details), but this is not needed at this time.

#### Step 6: Extract z-scores from all univariate GWAS (permuted and unpermuted)

#### Step 7: Run MOSTest on z-scores 

#### Step 8: Convert to csv
Convert the .mat output into csv using pvals2csv.py script

#### Step 9: Clump sumstats
Clump the results with sumstats.py

#### Nextflow pipeline for steps 3-8:
The nextflow pipeline "mostest_with_nextflow_v1.nf" can be used to complete step 3 to 8 in one go. Read, modify and submit the "implement_nextflow_mostest_v1.job" to accomplish this pipeline. 
