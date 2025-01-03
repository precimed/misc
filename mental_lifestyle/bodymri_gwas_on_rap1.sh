dx cat nt/project/norment/metabolimics/nmr_cauc_sep23_35_most_orig_clumpled.lead.snp2gene.csv | awk 'BEGIN { FS = "\t" } ; {print $8}' | sed '1d' > nmr_gene_sets.txt

#ls

#grep -F -f nmr_gene_sets.txt gene_sets_all.txt > nmr_gene_sets_ready.txt

#wc -l *txt

#head nmr_gene_sets_ready.txt
#project-G64PFY0JyZk3Gq1xBQYfk45Q:/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/
# Specify phenotype file and BGEN genotype file prefix
phenotype_all_sample="/zillur/bodymri/pheno/phenotype_bodymri_caucasian_all_sample_transformed_fin1.csv"
phenotype_only_male="/zillur/bodymri/pheno/phenotype_bodymri_caucasian_male_transformed_fin1.csv"
phenotype_only_female="/zillur/bodymri/pheno/phenotype_bodymri_caucasian_female_transformed_fin1.csv"
covariate_file1="/zillur/bodymri/pheno/covariates_age_sex_bmi_pc10_caucasian.csv"
covariate_file2="/zillur/bodymri/pheno/covariates_age_sex_pc10_caucasian.csv"
covariate_file3="/zillur/bodymri/pheno/covariates_age_nosex_bmi_pc10_caucasian.csv"
covariate_file4="/zillur/bodymri/pheno/covariates_age_nosex_pc10_caucasian.csv"
#ready_snp_list="project-G64PFY0JyZk3Gq1xBQYfk45Q:/norment/metabolimics/gene_sets_all.txt"
genotype_prefix="/ukb_c1-22_genome.filtered3_merged"
#head $phenotype_file
#project-G64PFY0JyZk3Gq1xBQYfk45Q:/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release/
bgen_prefix="/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release"
#ldl_custom="project-G64PFY0JyZk3Gq1xBQYfk45Q:/norment/metabolimics/ldl_custom_masks.txt"

#wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20231011.zip
#unzip plink2_linux_avx2_20231011.zip
#./plink2


data_field_bgen="ukb23159"
data_out1="/zillur/bodymri/gwas/qc/"
data_out2="/zillur/bodymri/gwas/step1/"
data_out3="/zillur/bodymri/gwas/step2/"
data_out4="/zillur/bodymri/gwas/burden/"

run_plink_all="plink2 --bfile ukb_c1-22_genome.filtered3_merged \
 --no-pheno --keep phenotype_bodymri_caucasian_all_sample_transformed_fin1.csv --autosome \
 --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 \
 --write-snplist --write-samples --no-id-header \
 --out qc_all_samples_pass --threads 32 --memory 16000"

run_plink_male="plink2 --bfile ukb_c1-22_genome.filtered3_merged \
 --no-pheno --keep phenotype_bodymri_caucasian_male_transformed_fin1.csv --autosome \
 --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 \
 --write-snplist --write-samples --no-id-header \
 --out qc_male_samples_pass --threads 32 --memory 16000"

run_plink_female="plink2 --bfile ukb_c1-22_genome.filtered3_merged \
 --no-pheno --keep phenotype_bodymri_caucasian_female_transformed_fin1.csv --autosome \
 --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 \
 --write-snplist --write-samples --no-id-header \
 --out qc_female_samples_pass --threads 32 --memory 16000"
 
#dx run swiss-army-knife -iin="${genotype_prefix}.bed" \
# -iin="${genotype_prefix}.bim" \
# -iin="${genotype_prefix}.fam" \
# -iin="${phenotype_all_sample}" \
# -icmd="${run_plink_all}" --tag="extract_good_snps" --instance-type "mem2_ssd1_v2_x32" \
# --destination="${project}:${data_out1}" --brief --yes

#dx run swiss-army-knife -iin="${genotype_prefix}.bed" \
# -iin="${genotype_prefix}.bim" \
# -iin="${genotype_prefix}.fam" \
# -iin="${phenotype_only_male}" \
# -icmd="${run_plink_male}" --tag="extract_good_snps" --instance-type "mem2_ssd1_v2_x32" \
# --destination="${project}:${data_out1}" --brief --yes

#dx run swiss-army-knife -iin="${genotype_prefix}.bed" \
# -iin="${genotype_prefix}.bim" \
# -iin="${genotype_prefix}.fam" \
# -iin="${phenotype_only_female}" \
# -icmd="${run_plink_female}" --tag="extract_good_snps" --instance-type "mem2_ssd1_v2_x32" \
# --destination="${project}:${data_out1}" --brief --yes

all_sample_snplist="/zillur/bodymri/gwas/qc/qc_all_samples_pass.snplist"
all_sample_id="/zillur/bodymri/gwas/qc/qc_all_samples_pass.id"
male_sample_snplist="/zillur/bodymri/gwas/qc/qc_male_samples_pass.snplist"
male_sample_id="/zillur/bodymri/qc/gwas/qc_male_samples_pass.id"
female_sample_snplist="/zillur/bodymri/gwas/qc/qc_female_samples_pass.snplist"
female_sample_id="/zillur/bodymri/gwas/qc/qc_female_samples_pass.id"

regenie_all_step1="regenie --step 1 \
 --lowmem --out step1_bodymri_all_with_bmi --bed ukb_c1-22_genome.filtered3_merged \
 --phenoFile phenotype_bodymri_caucasian_all_sample_transformed_fin1.csv \
 --covarFile covariates_age_sex_bmi_pc10_caucasian.csv \
 --extract qc_all_samples_pass.snplist \
 --bsize 1000 --loocv --gz --threads 96"

regenie_all_no_bmi_step1="regenie --step 1 \
 --lowmem --out step1_bodymri_all_no_bmi --bed ukb_c1-22_genome.filtered3_merged \
 --phenoFile phenotype_bodymri_caucasian_all_sample_transformed_fin1.csv \
 --covarFile covariates_age_sex_pc10_caucasian.csv \
 --extract qc_all_samples_pass.snplist \
 --bsize 1000 --loocv --gz --threads 96"

regenie_male_step1="regenie --step 1 \
 --lowmem --out step1_bodymri_male_with_bmi --bed ukb_c1-22_genome.filtered3_merged \
 --phenoFile phenotype_bodymri_caucasian_male_transformed_fin1.csv \
 --covarFile covariates_age_nosex_bmi_pc10_caucasian.csv \
 --extract qc_male_samples_pass.snplist \
 --bsize 1000 --loocv --gz --threads 96"

regenie_male_no_bmi_step1="regenie --step 1 \
 --lowmem --out step1_bodymri_male_no_bmi --bed ukb_c1-22_genome.filtered3_merged \
 --phenoFile phenotype_bodymri_caucasian_male_transformed_fin1.csv \
 --covarFile covariates_age_nosex_pc10_caucasian.csv \
 --extract qc_male_samples_pass.snplist \
 --bsize 1000 --loocv --gz --threads 96"

regenie_female_step1="regenie --step 1 \
 --lowmem --out step1_bodymri_female_with_bmi --bed ukb_c1-22_genome.filtered3_merged \
 --phenoFile phenotype_bodymri_caucasian_female_transformed_fin1.csv \
 --covarFile covariates_age_nosex_bmi_pc10_caucasian.csv \
 --extract qc_female_samples_pass.snplist \
 --bsize 1000 --loocv --gz --threads 96"

regenie_female_no_bmi_step1="regenie --step 1 \
 --lowmem --out step1_bodymri_female_no_bmi --bed ukb_c1-22_genome.filtered3_merged \
 --phenoFile phenotype_bodymri_caucasian_female_transformed_fin1.csv \
 --covarFile covariates_age_nosex_pc10_caucasian.csv \
 --extract qc_female_samples_pass.snplist \
 --bsize 1000 --loocv --gz --threads 96"

run_step2=(
"${run_regenie_gwas_all}"
"${run_regenie_gwas_all_no_bmi}"
"${run_regenie_gwas_male}"
"${run_regenie_gwas_male_no_bmi}"
"${run_regenie_gwas_female}"
"${run_regenie_gwas_female_no_bmi}"
)

run_burden=(
"${run_regenie_burden_all}"
"${run_regenie_burden_all_no_bmi}"
"${run_regenie_burden_male}"
"${run_regenie_burden_male_no_bmi}"
"${run_regenie_burden_female}"
"${run_regenie_burden_female_no_bmi}"
)
#for run_steps in "${run_step1[@]}"
#do
# dx run swiss-army-knife -iin="${genotype_prefix}.bed" \
#  -iin="${genotype_prefix}.bim" \
#  -iin="${genotype_prefix}.fam" \
#  -iin="${phenotype_all_sample}" -iin="${phenotype_only_male}" -iin="${phenotype_only_female}" \
#  -iin="${covariate_file1}" -iin="${covariate_file2}" -iin="${covariate_file3}" -iin="${covariate_file4}" \
#  -iin="${all_sample_snplist}" -iin="${male_sample_snplist}" -iin="${female_sample_snplist}" \
#  -icmd="${run_steps}" \
#  --tag="gwas_step1" --instance-type "mem1_hdd1_v2_x96" \
#  --destination="${project}:${data_out2}" --brief --yes
#done

pred_folder="/zillur/bodymri/gwas/step1"
path_to_500kwes_helper_files="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/helper_files"
ldl_custom="/zillur/norment/metabolimics/ldl_custom_masks.txt"
ready_snp_list="/zillur/norment/metabolimics/gene_sets_all.txt"

pred1=/zillur/bodymri/gwas/step1/step1_bodymri_all_with_bmi_pred.list
pred2=/zillur/bodymri/gwas/step1/step1_bodymri_all_no_bmi_pred.list
pred3=/zillur/bodymri/gwas/step1/step1_bodymri_male_with_bmi_pred.list
pred4=/zillur/bodymri/gwas/step1/step1_bodymri_male_no_bmi_pred.list
pred5=/zillur/bodymri/gwas/step1/step1_bodymri_female_with_bmi_pred.list
pred6=/zillur/bodymri/gwas/step1/step1_bodymri_female_no_bmi_pred.list

#dx run swiss-army-knife -icmd="cp /mnt/project/zillur/bodymri/gwas/step1/*loco.gz . ; zip loco.zip *loco.gz ; rm *loco.gz" --destination /zillur/bodymri/gwas/step1/ 


for i in {1..22}; do
    run_regenie_burden_all="regenie \
    --step 2 --pred step1_bodymri_all_with_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_all_sample_transformed_fin1.csv \
    --covarFile covariates_age_sex_bmi_pc10_caucasian.csv \
    --set-list ukb23158_500k_OQFE.sets.txt.gz \
    --anno-file ukb23158_500k_OQFE.annotations.txt.gz \
    --mask-def ldl_custom_masks.txt \
    --aaf-bins 0.005 --nauto 23 \
    --bsize 200 --extract-sets gene_sets_all.txt \
    --write-mask-snplist --write-mask \
    --vc-tests skato,acato-full \
    --out bodymri_all_gene_burden_chr${i}"
    
    run_regenie_burden_all_no_bmi="regenie \
    --step 2 --pred step1_bodymri_all_no_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_all_sample_transformed_fin1.csv \
    --covarFile covariates_age_sex_pc10_caucasian.csv \
    --set-list ukb23158_500k_OQFE.sets.txt.gz \
    --anno-file ukb23158_500k_OQFE.annotations.txt.gz \
    --mask-def ldl_custom_masks.txt \
    --aaf-bins 0.005 --nauto 23 \
    --bsize 200 --extract-sets gene_sets_all.txt \
    --write-mask-snplist --write-mask \
    --vc-tests skato,acato-full \
    --out bodymri_all_no_bmi_gene_burden_chr${i}"
    
    run_regenie_burden_male="regenie \
    --step 2 --pred step1_bodymri_male_with_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_male_transformed_fin1.csv \
    --covarFile covariates_age_nosex_bmi_pc10_caucasian.csv \
    --set-list ukb23158_500k_OQFE.sets.txt.gz \
    --anno-file ukb23158_500k_OQFE.annotations.txt.gz \
    --mask-def ldl_custom_masks.txt \
    --aaf-bins 0.005 --nauto 23 \
    --bsize 200 --extract-sets gene_sets_all.txt \
    --write-mask-snplist --write-mask \
    --vc-tests skato,acato-full \
    --out bodymri_male_gene_burden_chr${i}"
    
    run_regenie_burden_male_no_bmi="regenie \
    --step 2 --pred step1_bodymri_male_no_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_male_transformed_fin1.csv \
    --covarFile covariates_age_nosex_pc10_caucasian.csv \
    --set-list ukb23158_500k_OQFE.sets.txt.gz \
    --anno-file ukb23158_500k_OQFE.annotations.txt.gz \
    --mask-def ldl_custom_masks.txt \
    --aaf-bins 0.005 --nauto 23 \
    --bsize 200 --extract-sets gene_sets_all.txt \
    --write-mask-snplist --write-mask \
    --vc-tests skato,acato-full \
    --out bodymri_male_no_bmi_gene_burden_chr${i}"
    
    run_regenie_burden_female="regenie \
    --step 2 --pred step1_bodymri_female_with_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_female_transformed_fin1.csv \
    --covarFile covariates_age_nosex_bmi_pc10_caucasian.csv \
    --set-list ukb23158_500k_OQFE.sets.txt.gz \
    --anno-file ukb23158_500k_OQFE.annotations.txt.gz \
    --mask-def ldl_custom_masks.txt \
    --aaf-bins 0.005 --nauto 23 \
    --bsize 200 --extract-sets gene_sets_all.txt \
    --write-mask-snplist --write-mask \
    --vc-tests skato,acato-full \
    --out bodymri_female_gene_burden_chr${i}"
    
    run_regenie_burden_female_no_bmi="regenie \
    --step 2 --pred step1_bodymri_female_no_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_female_transformed_fin1.csv \
    --covarFile covariates_age_nosex_pc10_caucasian.csv \
    --set-list ukb23158_500k_OQFE.sets.txt.gz \
    --anno-file ukb23158_500k_OQFE.annotations.txt.gz \
    --mask-def ldl_custom_masks.txt \
    --aaf-bins 0.005 --nauto 23 \
    --bsize 200 --extract-sets gene_sets_all.txt \
    --write-mask-snplist --write-mask \
    --vc-tests skato,acato-full \
    --out bodymri_female_no_bmi_gene_burden_chr${i}"
    
    run_regenie_gwas_all="regenie \
    --step 2 --pred step1_bodymri_all_with_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_all_sample_transformed_fin1.csv \
    --covarFile covariates_age_sex_bmi_pc10_caucasian.csv \
    --bsize 400 --minINFO 0.8 --minMAC 20 \
    --out bodymri_all_step2_chr${i}"
    
    run_regenie_gwas_all_no_bmi="regenie \
    --step 2 --pred step1_bodymri_all_no_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_all_sample_transformed_fin1.csv \
    --covarFile covariates_age_sex_pc10_caucasian.csv \
    --bsize 400 --minINFO 0.8 --minMAC 20 \
    --out bodymri_all_no_bmi_step2_chr${i}"
    
    run_regenie_gwas_male="regenie \
    --step 2 --pred step1_bodymri_male_with_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_male_transformed_fin1.csv \
    --covarFile covariates_age_nosex_bmi_pc10_caucasian.csv \
    --bsize 400 --minINFO 0.8 --minMAC 20 \
    --out bodymri_male_step2_chr${i}"
    
    run_regenie_gwas_male_no_bmi="regenie \
    --step 2 --pred step1_bodymri_male_no_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_male_transformed_fin1.csv \
    --covarFile covariates_age_nosex_pc10_caucasian.csv \
    --bsize 400 --minINFO 0.8 --minMAC 20 \
    --out bodymri_male_no_bmi_step2_chr${i}"
    
    run_regenie_gwas_female="regenie \
    --step 2 --pred step1_bodymri_female_with_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_female_transformed_fin1.csv \
    --covarFile covariates_age_nosex_bmi_pc10_caucasian.csv \
    --bsize 400 --minINFO 0.8 --minMAC 20 \
    --out bodymri_female_step2_chr${i}"
    
    run_regenie_gwas_female_no_bmi="regenie \
    --step 2 --pred step1_bodymri_female_no_bmi_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile phenotype_bodymri_caucasian_female_transformed_fin1.csv \
    --covarFile covariates_age_nosex_pc10_caucasian.csv \
    --bsize 400 --minINFO 0.8 --minMAC 20 \
    --out bodymri_female_no_bmi_step2_chr${i}"
   
    dx run swiss-army-knife -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.bgen" \
     -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.sample" \
     -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.bgen.bgi" \
     -iin="${ldl_custom}" \
     -iin="${ready_snp_list}" \
     -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.sets.txt.gz" \
     -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.annotations.txt.gz" \
     -iin="${phenotype_all_sample}" -iin="${phenotype_only_male}" -iin="${phenotype_only_female}" \
     -iin="${covariate_file1}" -iin="${covariate_file2}" -iin="${covariate_file3}" -iin="${covariate_file4}" \
     -iin="${pred1}" -iin="${pred2}" -iin="${pred3}" -iin="${pred4}" -iin="${pred5}" -iin="${pred6}" \
     -iin="${pred_folder}/loco.zip" \
     -icmd="unzip loco.zip ;  ${run_regenie_gwas_all} ; ${run_regenie_gwas_all_no_bmi} ; \
       ${run_regenie_gwas_male} ; ${run_regenie_gwas_male_no_bmi} ; \
       ${run_regenie_gwas_female} ; ${run_regenie_gwas_female_no_bmi}" \
     --tag="bodymri_step2" --instance-type "mem2_ssd1_v2_x96" \
     --destination="${project}:${data_out3}" --brief --yes
     
     dx run swiss-army-knife -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.bgen" \
      -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.sample" \
      -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.bgen.bgi" \
      -iin="${ldl_custom}" \
      -iin="${ready_snp_list}" \
      -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.sets.txt.gz" \
      -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.annotations.txt.gz" \
      -iin="${phenotype_all_sample}" -iin="${phenotype_only_male}" -iin="${phenotype_only_female}" \
      -iin="${covariate_file1}" -iin="${covariate_file2}" -iin="${covariate_file3}" -iin="${covariate_file4}" \
      -iin="${pred1}" -iin="${pred2}" -iin="${pred3}" -iin="${pred4}" -iin="${pred5}" -iin="${pred6}" \
      -iin="${pred_folder}/loco.zip" \
      -icmd="unzip loco.zip ; ${run_regenie_burden_all} ; ${run_regenie_burden_all_no_bmi} ; \
       ${run_regenie_burden_male} ; ${run_regenie_burden_male_no_bmi} ; \
       ${run_regenie_burden_female} ; ${run_regenie_burden_female_no_bmi}" \
      --tag="bodymri_burden" --instance-type "mem2_ssd1_v2_x96" \
      --destination="${project}:${data_out4}" --brief --yes
done  



