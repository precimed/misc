path_to_500kwes_helper_files="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/helper_files"
#dx ls "${path_to_500kwes_helper_files}"

#dx zcat "${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.sets.txt.gz" | awk '{print $1}' > gene_sets_all.txt

#dx ls /mnt/project/norment/metabolimics

#dx cat nt/project/norment/metabolimics/nmr_cauc_sep23_35_most_orig_clumpled.lead.snp2gene.csv | awk 'BEGIN { FS = "\t" } ; {print $8}' | sed '1d' > nmr_gene_sets.txt

#ls

#grep -F -f nmr_gene_sets.txt gene_sets_all.txt > nmr_gene_sets_ready.txt

#wc -l *txt

#head nmr_gene_sets_ready.txt
#project-G64PFY0JyZk3Gq1xBQYfk45Q:/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/
# Specify phenotype file and BGEN genotype file prefix
phenotype_file="/norment/metabolimics/ukbnmr/burden/lipid_pheno_dec_with_labels_tech_var_norelative_with_NA.csv"
covariate_file="/norment/metabolimics/ukbnmr/burden/covariate_sep23_age_sex_20pc.csv"
ready_snp_list="project-G64PFY0JyZk3Gq1xBQYfk45Q:/norment/metabolimics/gene_sets_all.txt"
genotype_prefix="/ukb_c1-22_genome.filtered3_merged"
#head $phenotype_file
#project-G64PFY0JyZk3Gq1xBQYfk45Q:/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release/
bgen_prefix="project-G64PFY0JyZk3Gq1xBQYfk45Q:/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release"
ldl_custom="project-G64PFY0JyZk3Gq1xBQYfk45Q:/norment/metabolimics/ldl_custom_masks.txt"

#wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20231011.zip
#unzip plink2_linux_avx2_20231011.zip
#./plink2

data_field_geno="ukb23158"
data_field_bgen="ukb23159"
data_out1="/norment/metabolimics/ukbnmr/burden/results/"

#data_out1="/norment/metabolimics/burden/good_snps/"

#run_plink_wes="plink2 --bfile ukb_c1-22_genome.filtered3_merged \
# --no-pheno --keep lipid_pheno_dec_with_labels_tech_var_norelative_with_NA.csv --autosome \
# --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 \
# --write-snplist --write-samples --no-id-header \
# --out WES_c1_22_snps_qc_pass --threads 32 --memory 16000"
 
#dx run swiss-army-knife -iin="${genotype_prefix}.bed" \
# -iin="${genotype_prefix}.bim" \
# -iin="${genotype_prefix}.fam" \
# -iin="${phenotype_file}" \
# -icmd="${run_plink_wes}" --tag="extract_good_snps" --instance-type "mem2_ssd1_v2_x32" \
# --destination="${project}:${data_out1}" --brief --yes

#run_regenie_step1="regenie --step 1 \
# --lowmem --out step1_metabolomics --bed ukb_c1-22_genome.filtered3_merged \
# --phenoFile lipid_pheno_dec_with_labels_tech_var_norelative_with_NA.csv \
# --covarFile covariate_sep23_age_sex_20pc.csv \
# --extract WES_c1_22_snps_qc_pass.snplist \
# --bsize 1000 --loocv --gz --threads 96"

#dx run swiss-army-knife -iin="${genotype_prefix}.bed" \
# -iin="${genotype_prefix}.bim" \
# -iin="${genotype_prefix}.fam" \
# -iin="${phenotype_file}" \
# -iin="${covariate_file}" \
# -iin="${data_out1}WES_c1_22_snps_qc_pass.snplist" \
# -icmd="${run_regenie_step1}" --tag="Step1" --instance-type "mem1_hdd1_v2_x96" \
# --destination="${project}:${data_out1}" --brief --yes


data_field_bgen="ukb23159"
pred_folder="/norment/metabolimics/ukbnmr/burden/snplist"

for i in {1..22}; do
    run_regenie_burden="regenie \
    --step 2 --pred step1_metabolomics_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile lipid_pheno_dec_with_labels_tech_var_norelative_with_NA.csv \
    --covarFile covariate_sep23_age_sex_20pc.csv \
    --set-list ukb23158_500k_OQFE.sets.txt.gz \
    --anno-file ukb23158_500k_OQFE.annotations.txt.gz \
    --mask-def ldl_custom_masks.txt \
    --aaf-bins 0.005 --nauto 23 \
    --bsize 200 --extract-sets gene_sets_all.txt \
    --write-mask-snplist --write-mask \
    --vc-tests skato,acato-full \
    --out ukbnmr_all_gene_burden_chr${i}"
    
 dx run swiss-army-knife -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.bgen" \
  -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.sample" \
  -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.bgen.bgi" \
  -iin="${phenotype_file}" \
  -iin="${covariate_file}" \
  -iin="${ldl_custom}" \
  -iin="${ready_snp_list}" \
  -iin="${pred_folder}/step1_metabolomics_pred.list" \
  -iin="${pred_folder}/loco.zip" \
  -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.sets.txt.gz" \
  -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.annotations.txt.gz" \
  -icmd="unzip loco.zip ; ${run_regenie_burden}" --tag="burden_ukbnmr" --instance-type "mem2_ssd1_v2_x96" \
  --destination="${project}:${data_out1}" --brief --yes
done
