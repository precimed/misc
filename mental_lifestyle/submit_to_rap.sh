pheno_file=/zillur/bip_bmi_smd/pheno/proteomic_pheno_caucasian_imputed_filtered_fin2.csv
dx upload run_caret_on_rap.R
run_caret_step1="Rscript run_caret_on_rap.R"
data_out1=/zillur/bip_bmi_smd/caret/
dx run swiss-army-knife -iin="run_caret_on_rap.R" \
 -iin="${pheno_file}" \
 -icmd="${run_caret_step1}" --tag="Caret" --instance-type "mem3_ssd1_v2_x96" \
 --destination="${project}:${data_out1}" --brief --yes

