//Lets declare inputs via channels

bgen_files=Channel.fromFilePairs("/cluster/projects/p33/s3-api/ukblake/genetics/imp/ukb_imp_chr*_v3.bgen",size:1).map {group_key,file_list -> tuple(group_key.replaceFirst(/^*ukb_imp_/,""),file_list.first())}
sample_files=Channel.fromPath("/cluster/projects/p33/s3-api/ukblake/genetics/imp/ukb27412_imp_chr1_v3_s487395.sample")
Channel.fromPath("${params.pheno}").into{pheno1;pheno2;pheno3}
Channel.fromPath("${params.cov}").into{cov1;cov2;cov3}
merging_list=Channel.fromPath("/cluster/projects/p33/users/mohammadzr/metabolomics/liver/pheno/merge_list_for_relatedness.txt")
unrelated_ind=Channel.from("unrelated_geno")

process merge_1 {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
module 'plink/1.90b6.2' 
executor='slurm'
queueSize='22'
jobName='merging'
cpus='16'

clusterOptions  "-A p33 -t 192:00:00 --mem-per-cpu 16000M"

publishDir params.out, mode:'copy'

input:
file merge_in from merging_list
file pheno from pheno1


output:
set path("*.bed"), path("*.bim"), path("*.fam") into merged

script:
outPrefix = 'merged_geno'

"""
plink --merge-list ${merge_in} --threads 16 --keep ${pheno} --make-bed --out ${outPrefix}
"""
}

process unrelate_1 {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
module 'plink2/2.00a2LM'
executor='slurm'
queueSize='22'
jobName='unrel'
cpus='16'
clusterOptions  "-A p33 -t 242:00:00 --mem-per-cpu 16000M"

publishDir params.out, mode:'copy'

input:
file merge_in2 from merged

output:
path("*.fam") into unrelated1,unrelated2

script:
outPrefix = 'unrelated_individual'
inprefix = 'merged_geno'
"""
plink2 --bfile ${inprefix} --king-cutoff 0.05 --make-just-fam --out ${outPrefix} --threads 16
"""
}

Channel.fromList(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22']).into {chrom1;chrom2}

bgen_files.combine(sample_files).combine(unrelated1).groupTuple().set {bgen1}
chrom1.join(bgen1).into{bgen2;bgen3;bgen4}

process extract_1 {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
module 'plink2/2.00a2LM'
executor='slurm'
queueSize='22'
jobName='extraction'
cpus='16'
clusterOptions  "-A p33 -t 192:00:00 --mem-per-cpu 16000M"

publishDir params.out, mode:'copy'

input:
tuple val(chr),path(bgen),path(sample),path(fam) from bgen2

output:
tuple val("${chr}"),path("${chr}*.*") optional true into extracted_imp1,extracted_imp2,extracted_imp3

script:

"""
plink2 --bgen ${bgen} ref-first --sample ${sample} --keep ${fam} --mac 20 --make-bed --out "${chr}_unrelated_imp" --threads 16
"""
}

process make_chunks {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
executor='slurm'
queueSize='22'
jobName='chunk_strt'
cpus='2'
clusterOptions  "-A p33 -t 24:00:00 --mem-per-cpu 8000M"

publishDir params.out, mode:'copy'

input:
tuple val(chr),path(bim) from extracted_imp1

output:
tuple val("${chr}"),path("${chr}*.*") optional true into chunks1

script:
"""
python /cluster/projects/p33/users/mohammadzr/metabolomics/scripts/make_chunks_by_snps.py ${chr}_unrelated_imp.bim
"""
}
chunks2=chunks1.transpose()
chunks2.combine(extracted_imp2,by:0).set {chunks3}

process chunks_to_plink {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
module 'plink/1.90b6.2'
executor='slurm'
queueSize='3456'
jobName='chnk2plink'
cpus='8'
clusterOptions  "-A p33 -t 24:00:00 --mem-per-cpu 8000M"

publishDir params.out, mode:'copy'

input:
tuple val(chr),path(chunk),path(bfiles) from chunks3

output:
tuple val("${chr}"),path("${chr}*.*") optional true into chunks4,chunks5

script:
def (bed,bim,fam,log) = bfiles
"""
plink --bed ${bed} --bim ${bim} --fam ${fam} --extract ${chunk} --make-bed --out "${chunk}" --threads 8
"""
}
chunks4.map{chr,bfiles -> tuple(chr,bfiles)}.combine(cov1).combine(pheno2).into{gwas_input1;gwas_input2}

process gwas1 {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
module 'plink2/2.00a2LM'
executor='slurm'
queueSize='3456'
jobName='gwas1'
cpus='16'
clusterOptions  "-A p33 -t 192:00:00 --mem-per-cpu 16000M"

publishDir params.out, mode:'copy'

input:
tuple val(chr),path(bfiles),path(cov),path(pheno) from gwas_input1

output:
tuple val("${chr}"),path("${chr}*.*") optional true into first_gwas

script:
def (beds,bim,fam,log)=bfiles
bed=beds.name.split('\\.')[0]
"""
plink2 --bfile ${bed}.csv --glm omit-ref hide-covar --covar ${cov} --covar-variance-standardize --pheno ${pheno} --out ${bed}_glm --threads 16
"""
}

process permutation {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
executor='slurm'
queueSize='3456'
jobName='permute'
cpus='2'
clusterOptions  "-A p33 -t 192:00:00 --mem-per-cpu 8000M"

publishDir params.out, mode:'copy'

input:
tuple val(chr),path(bfiles) from chunks5

output:
tuple val("${chr}"),path("${chr}*.*") optional true into permuted

script:
def (beds,bim,fam,log)=bfiles
bed=beds.name.split('\\.')[0]
"""
python /cluster/projects/p33/users/alexeas/mostest_ukb/scripts/permute_bed.py --bfile ${bed}.csv --out ${bed}_permuted
"""
}

permuted.combine(cov2).combine(pheno3).into{permuted_input1;permuted_input2}

process gwas2 {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
module 'plink2/2.00a2LM'
executor='slurm'
queueSize='3456'
jobName='gwas2'
cpus='16' 
clusterOptions  "-A p33 -t 192:00:00 --mem-per-cpu 16000M"

publishDir params.out, mode:'copy'

input:
tuple val(chr),path(bfiles),path(cov),path(pheno) from permuted_input1

output:
tuple val("${chr}"),path("${chr}*.*") optional true into second_gwas

script:
def (beds,bim,fam,log)=bfiles
bed=beds.name.split('\\.')[0]
"""
plink2 --bfile ${bed} --glm omit-ref hide-covar --covar ${cov} --covar-variance-standardize --pheno ${pheno} --out ${bed} --threads 16
"""
}
first_gwas.transpose().map{chr,trait -> tuple(trait.name.split('\\.')[1],trait,chr)}.groupTuple().into{merge_in1;merge_in2}

process merge_chunks {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
executor='slurm'
queueSize='200'
jobName='merge_ch1'
cpus='2'
clusterOptions  "-A p33 -t 24:00:00 --mem-per-cpu 8000M"

publishDir params.out, mode:'copy'

input:
tuple val(trait),path(sumstats),val(chr) from merge_in1

output:
tuple val("${trait}"),path("${trait}*.*") optional true into merge_out1

script:
"""
if [[ $trait != 'log' ]]; then
sed -e '1,\${/^\\#/d' -e '}' ${sumstats} | sort -k1,1n -k2,2n > ${trait}_original.csv
fi
"""
}
merge_out1.map{trt,smt -> tuple(smt.name.split('\\.')[1],smt,trt)}.groupTuple().set{z_in1}
z_in1.map{csv,smt,trt -> tuple(csv,smt.size(),smt,trt)}.set{z_in3}

second_gwas.transpose().map{chr,trait -> tuple(trait.name.split('\\.')[1],trait,chr)}.groupTuple().into{merge_in3;merge_in4}

process merge_chunks2 {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
executor='slurm'
queueSize='200'
jobName='merge_ch2'
cpus='2'
clusterOptions  "-A p33 -t 24:00:00 --mem-per-cpu 8000M"

publishDir params.out, mode:'copy'

input:
tuple val(trait),path(sumstats),val(chr) from merge_in3

output:
tuple val("${trait}"),path("${trait}*.*") optional true into merge_out2,merge_out3

script:
"""
if [[ $trait != 'log' ]]; then
sed -e '1,\${/^\\#/d' -e '}' ${sumstats} | sort -k1,1n -k2,2n > ${trait}_permuted.csv
fi
"""
}
merge_out2.map{trt,smt -> tuple(smt.name.split('\\.')[1],smt,trt)}.groupTuple().set{z_in2}
z_in2.map{csv,smt,trt -> tuple(csv,smt.size(),smt,trt)}.set{z_in4}

process merge_zscore1 {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
executor='slurm'
queueSize='200'
jobName='merge_z'
cpus='2'
clusterOptions  "-A p33 -t 24:00:00 --mem-per-cpu 8000M"

publishDir params.out, mode:'copy'

input:
tuple val(csv1),val(sz1),path(sumstats1),val(trait1) from z_in3
output:
tuple path("glm_original_gwas_tstat.csv"),path("glm_original_gwas_zscore.csv") optional true into z_out1,z_m1

script:
"""
ARG=\$(awk 'BEGIN{ORS=""; N=11; print(N); for(i=1; i<${sz1}; i++) {N=N+13; print(","N)}}')
paste ${sumstats1} | cut -f \$ARG  > glm_original_gwas_tstat.csv
echo ${trait1} | sed -e 's/,/\\t/g;s/\\[\\|\\]//g' | cat - glm_original_gwas_tstat.csv > glm_original_gwas_zscore.csv
"""
}
process merge_zscore2 {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
executor='slurm'
queueSize='200'
jobName='merge_z'
cpus='2'
clusterOptions  "-A p33 -t 24:00:00 --mem-per-cpu 8000M"

publishDir params.out, mode:'copy'

input:
tuple val(csv2),val(sz2),path(sumstats2),val(trait2) from z_in4
output:
tuple path("glm_permuted_gwas_tstat.csv"),path("glm_permuted_gwas_zscore.csv") optional true into z_out2

script:
"""
ARG=\$(awk 'BEGIN{ORS=""; N=11; print(N); for(i=1; i<${sz2}; i++) {N=N+13; print(","N)}}')
paste ${sumstats2} | cut -f \$ARG  > glm_permuted_gwas_tstat.csv
echo ${trait2} | sed -e 's/,/\\t/g;s/\\[\\|\\]//g' | cat - glm_permuted_gwas_tstat.csv > glm_permuted_gwas_zscore.csv
"""
}

process prep {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
module 'R/4.2.0-foss-2021b'
executor='slurm'
queueSize='2'
jobName='mostest'
cpus='8'
clusterOptions  "-A p33 -t 24:00:00 --mem-per-cpu 8000M"

publishDir params.out, mode:'copy'

input:
tuple path(tstat1),path(z1) from z_out1
tuple path(tstat2),path(z2) from z_out2
output:
tuple path("${z1}_ordered_cols.csv"),path("${z2}_ordered_cols.csv") optional true into z_m11,z_m12,z_m13

script:
"""
Rscript /cluster/projects/p33/users/mohammadzr/metabolomics/scripts/order_zscore_cols.R ${z1} ${z2}
"""
}
//z_m11.view()
z_m11.map{it[0]}.splitCsv(sep:'\t',limit:1).map{(1..it.size()).collect()}.set{m1}
m1.transpose().flatten().combine(z_m12).into{m2;m3;m4}
z_m13.map{it[0]}.splitCsv(sep:'\t',limit:1).map{(it.size())}.set{m14}
channel.of(params.project).flatten().set{pro1}
m4.combine(pro1).combine(m14).set{m5}

process mostest {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
module 'MATLAB/2020b'
executor='slurm'
queueSize='2'
jobName='mostest'
cpus='16'
clusterOptions  "-A p33 -t 2:00:00 --mem-per-cpu 16000M"

publishDir params.out, mode:'copy'

input:
tuple val(a1),path(z1),path(z2),val(project1),val(lim) from m5

output:
tuple val("${a1}"),path("${my_mat}_${a1}.*") optional true into most

script:
"""
if [[ $a1 < $lim ]]; then
matlab -nodisplay -nosplash -batch "zmat_orig_file='${z1}';zmat_perm_file='${z2}'; num_eigval_to_keep=${a1}; out='${project1}_${a1}'; run /cluster/projects/p33/users/mohammadzr/metabolomics/scripts/mostest_code/mostest_mental; exit"
fi
"""
}
merge_out3.map{it[1]}.flatten().first().set{in_n1}
extracted_imp3.map{it[1]}.map{it[1]}.flatten().toSortedList().set{in_bim1}

process get_n {
beforeScript 'source /cluster/projects/p33/users/mohammadzr/envs/nextf/bin/activate'
executor='slurm'
queueSize='2'
jobName='mostest'
cpus='2'
clusterOptions  "-A p33 -t 24:00:00 --mem-per-cpu 8000M"

publishDir params.out, mode:'copy'

input:
file(n1) from in_n1
file(b1) from in_bim1
output:
tuple path("variant_sample.csv"),path("bim_for_convert.bim") optional true into bim_and_n

script:
"""
awk -v OFS='\t' '{print \$3,\$8}' ${n1} > variant_sample.csv
cat ${b1} | sort -k1,1n -k2,4n > bim_for_convert.bim
"""
}
