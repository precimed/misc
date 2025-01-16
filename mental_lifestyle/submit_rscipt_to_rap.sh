pheno_file=/zillur/bip_bmi_smd/pheno/proteomic_pheno_caucasian_imputed_filtered_fin2.csv
dx upload run_caret_on_rap3.R --path=/zillur/bip_bmi_smd/caret/
run_caret_step1="Rscript run_caret_on_rap3.R proteomic_pheno_caucasian_imputed_filtered_fin2.csv"
my_script=/zillur/bip_bmi_smd/caret/run_caret_on_rap3.R
data_out1=/zillur/bip_bmi_smd/caret/
dx run swiss-army-knife -iin="${pheno_file}"\
 -iin="${my_script}" \
 -icmd="${run_caret_step1}" --tag="Caret" --instance-type "mem3_ssd1_v2_x96" \
 --destination="${project}:${data_out1}" --brief --yes

