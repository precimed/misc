library(ieugwasr)
library(foreach)
library(doParallel)
library(readr)  # Make sure you load the readr library if needed
library(TwoSampleMR)
library(MRPRESSO)
library(tidyverse)
library(parallel)
args = commandArgs(trailingOnly=TRUE)

ld_clump_local=function (dat, clump_kb, clump_r2, clump_p, bfile, plink_bin,memory) 
{
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", 
                  "sh")
  fn <- tempfile()
  write.table(data.frame(SNP = dat[["rsid"]], P = dat[["pval"]]), 
              file = fn, row.names = F, col.names = T, quote = F)
  fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --clump ", shQuote(fn, 
                                                                    type = shell), " --clump-p1 ", clump_p, " --clump-r2 ", 
                 clump_r2, " --clump-kb ", clump_kb, " --memory ", memory, " --out ", shQuote(fn, 
                                                                                              type = shell))
  system(fun2)
  res <- read.table(paste(fn, ".clumped", sep = ""), header = T)
  unlink(paste(fn, "*", sep = ""))
  y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
  if (nrow(y) > 0) {
    message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), 
            " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}

exposure_metabolites=function(sumstats){
  exposure=list.files(sumstats,full.names = T)
  exposure1 <- mclapply(exposure,function(file1){
    tr1=strsplit(file1,"/")[[1]][[11]]
    tr2=gsub('_combined_original.csv','',tr1)
    start_time <- Sys.time()
    system(paste("echo Exposure: ",tr2," Started at: " , start_time,"\n"))
    flush.console()
    data = read.table(file1, sep="\t", header = T,
                      stringsAsFactors = F,comment.char = "")
    colnames(data)[1]='CHROM'
    exposure_gws1<-filter(data, P < 5e-08)
    exposure_gws_format1<-format_data(exposure_gws1, snp_col = "ID",
                                      beta_col = "BETA", effect_allele_col = "A1",se_col = "SE",
                                      other_allele_col = "REF",pval_col = "P",
                                      samplesize_col = "OBS_CT",chr_col = "CHROM", pos_col = "POS")
    exposure_gws_format1$exposure = tr2
    exposure_gws_format_clump1<-ld_clump_local(dplyr::tibble(rsid=exposure_gws_format1$SNP, 
                                                             pval=exposure_gws_format1$pval.exposure,
                                                             id=exposure_gws_format1$exposure),
                                               plink_bin = "/cluster/projects/p33/users/mohammadzr/metabolomics/1k_data/plink",
                                               bfile = "/cluster/projects/p33/users/mohammadzr/metabolomics/1k_data/EUR",
                                               clump_p = 1,clump_r2 = 0.001,clump_kb = 10000,memory = 1800)
    exposure_gws_format_clump1=as.data.frame(exposure_gws_format_clump1)
    exposure_gws_format_clump2=exposure_gws_format1[exposure_gws_format1$SNP %in% exposure_gws_format_clump1$rsid,]
    exposure_gws_format_clump2$pval.exposure=exposure_gws_format_clump1$pval[match(exposure_gws_format_clump1$rsid,
                                                                                   exposure_gws_format_clump2$SNP)]
    
    lv1=list(tr2,exposure_gws_format_clump2)
    end_time <- Sys.time()
    system(paste("echo Exposure: ",tr2," Took: ", end_time - start_time,"\n"))
    flush.console()
    return(lv1)},mc.cores = 128)
  #print(exposure1)
  return(exposure1)
}
#sumstats_metabolite='/cluster/projects/p33/users/mohammadzr/metabolomics/lipids/sep13/test'
sumstats_metabolite= args[1]

exposure_smds=function(sumstats){
  exposure=list.files(sumstats,full.names = T)
  exposure1 <- mclapply(exposure,function(file1){
    tr1=strsplit(strsplit(file1,'/')[[1]][9],'_')[[1]][2]
    start_time <- Sys.time()
    system(paste("echo Exposure: ",tr1," Started at: ", start_time,"\n"))
    flush.console()
    exposure1 <- read.table(file1, sep="\t", header = T,stringsAsFactors = F)
    exposure_gws1<-filter(exposure1, PVAL < 5e-08)
    if (tr1 == "MDD"){
      exposure_gws_format1<-format_data(exposure_gws1,snp_col = "SNP",beta_col = "BETA",
                                        effect_allele_col = "A1",se_col = "SE",other_allele_col = "A2",
                                        pval_col = "PVAL",samplesize_col = "N",chr_col = "CHR",
                                        pos_col = "BP")
      cat("BETA as beta_col :", tr1)
    } else {exposure_gws_format1<-format_data(exposure_gws1,snp_col = "SNP",beta_col = "OR",
                                              effect_allele_col = "A1",se_col = "SE",
                                              other_allele_col = "A2",pval_col = "PVAL",
                                              samplesize_col = "N",chr_col = "CHR", pos_col = "BP")
    cat("OR as beta_col: ", tr1)
    }
    exposure_gws_format1$exposure = tr1
    exposure_gws_format_clump1<-ld_clump_local(dplyr::tibble(rsid=exposure_gws_format1$SNP, 
                                                             pval=exposure_gws_format1$pval.exposure,
                                                             id=exposure_gws_format1$exposure),
                                               plink_bin = "/cluster/projects/p33/users/mohammadzr/metabolomics/1k_data/plink",
                                               bfile = "/cluster/projects/p33/users/mohammadzr/metabolomics/1k_data/EUR",
                                               clump_p = 1,clump_r2 = 0.001,clump_kb = 10000,memory = 1800)
    
    if (tr1 == "BIP" | tr1 == "SCZ"){
      exposure_gws_format_clump1=as.data.frame(exposure_gws_format_clump1)
      exposure_gws_format_clump2=exposure_gws_format1[exposure_gws_format1$SNP %in% exposure_gws_format_clump1$rsid,]
      exposure_gws_format_clump2$pval.exposure=exposure_gws_format_clump1$pval[match(exposure_gws_format_clump1$rsid,
                                                                                     exposure_gws_format_clump2$SNP)]
      exposure_gws_format_clump2$beta.exposure=log(exposure_gws_format_clump2$beta.exposure)
      cat("Exposure log: ", tr1)
    } else {exposure_gws_format_clump1=as.data.frame(exposure_gws_format_clump1)
    exposure_gws_format_clump2=exposure_gws_format1[exposure_gws_format1$SNP %in% exposure_gws_format_clump1$rsid,]
    exposure_gws_format_clump2$pval.exposure=exposure_gws_format_clump1$pval[match(exposure_gws_format_clump1$rsid,
                                                                                   exposure_gws_format_clump2$SNP)]
    
    cat("No log exposure: ",tr1)
    }
    lv1=list(tr1,exposure_gws_format_clump2)
    end_time <- Sys.time()
    system(paste("echo Exposure: ",tr1," Took: ", end_time - start_time,"\n"))
    flush.console()
    return(lv1)},mc.cores = 128)
  return(exposure1)
}
#sumstats_smd='/cluster/projects/p33/users/mohammadzr/metabolomics/smd_sumstats'
sumstats_smd=args[2]

smd_outcomes <- function(sumstats, exposures) {
  smd <- list.files(sumstats, full.names = TRUE)
  
  my_out <- lapply(exposures, function(exposure) {
    exposure_name <- exposure[[1]]
    exposure_data <- exposure[[2]]
    smd_outcome <- list()
    
    outcome_list <- mclapply(smd, function(file1) {
      trait_name <- strsplit(strsplit(file1, '/')[[1]][9], '_')[[1]][2]
      start_time <- Sys.time()
      system(paste("echo Exposure: ", exposure_name, " Outcome: ", trait_name,
                   " Started at: ", start_time,"\n"))
      flush.console()
      
      if (trait_name == "MDD") {
        outcome_format1 <- read_outcome_data(
          snps = exposure_data$SNP,filename = file1, sep = "\t",snp_col = "SNP",
          beta_col = "BETA", se_col = "SE",effect_allele_col = "A1",
          other_allele_col = "A2",pval_col = "PVAL",samplesize_col = "N")
      } else {
        outcome_format1 <- read_outcome_data(
          snps = exposure_data$SNP,filename = file1,sep = "\t",snp_col = "SNP",
          beta_col = "OR",se_col = "SE",effect_allele_col = "A1",
          other_allele_col = "A2",pval_col = "PVAL",samplesize_col = "N")
        outcome_format1$beta.outcome <- log(outcome_format1$beta.outcome)
      }
      outcome_format1$outcome <- trait_name
      H_data1 <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_format1)
      print(dim(H_data1))
      end_time <- Sys.time()
      elapsed_time <- end_time - start_time
      system(paste("echo Exposure: ", exposure_name, " Outcome: ", trait_name, " End at: ", end_time,
                   "and total time took: ", elapsed_time,"\n"))
      flush.console()
      return(H_data1)
    },mc.cores = 128)
    
    smd_outcome[[exposure_name]] <- outcome_list
    return(smd_outcome)
  })
  
  return(my_out)
}

metabolite_outcomes <- function(sumstats, exposures) {
  metabolite <- list.files(sumstats, full.names = TRUE)
  
  my_out <- lapply(exposures, function(exposure) {
    exposure_name <- exposure[[1]]
    exposure_data <- exposure[[2]]
    metabolite_outcome <- list()
    
    outcome_list <- mclapply(metabolite, function(file1) {
      tr1 <- strsplit(file1, "/")[[1]][[11]]
      tr2 <- gsub('_combined_original.csv', '', tr1)
      start_time <- Sys.time()
      system(paste("echo Exposure: ", exposure_name, " Outcome: ", tr2, " Started at: ", start_time,"\n"))
      flush.console()
      outcome_format1 <- read_outcome_data(
        snps = exposure_data$SNP,filename = file1,sep = "\t",snp_col = "ID",
        beta_col = "BETA",se_col = "SE",effect_allele_col = "A1",
        other_allele_col = "REF",pval_col = "P",samplesize_col = "OBS_CT")
      
      outcome_format1$outcome <- tr2
      H_data1 <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_format1)
      print(dim(H_data1))
      end_time <- Sys.time()
      elapsed_time <- end_time - start_time
      system(paste("echo Exposure: ", exposure_name, " Outcome: ", tr2, " End at: ", end_time,
                   "and total time took: ", elapsed_time,"\n"))
      flush.console()
      return(H_data1)
    },mc.cores = 128)
    
    metabolite_outcome[[exposure_name]] <- outcome_list
    return(metabolite_outcome)
  })
  
  return(my_out)
}


# Call the metabolite_outcomes function
exposure_metabolite=exposure_metabolites(sumstats_metabolite)
exposure_smd=exposure_smds(sumstats_smd)
metabolite_outcome <- metabolite_outcomes(sumstats_metabolite, exposure_smd)
mental_outcome <- smd_outcomes(sumstats_smd, exposure_metabolite)


mr_my_methods=c("mr_ivw", "mr_weighted_median", "mr_weighted_mode",
                "mr_egger_regression_bootstrap")

# Define a custom function to apply to data frames
my_mr <- function(df) {
  # Perform your operations on the data frame here
  # For example, let's add 1 to the 'beta.exposure' column
  df1 <- mr(df,method_list = mr_my_methods)
  return(df1)
}

mr_metabolite_outcome <- lapply(metabolite_outcome, function(sublist) {
  processed_list <- lapply(sublist, function(df) {
    processed_df <- mclapply(df, my_mr, mc.cores = 128)
    return(processed_df)
  })
  return(processed_list)
})

mr_smd_outcome <- lapply(mental_outcome, function(sublist) {
  processed_list <- lapply(sublist, function(df) {
    processed_df <- mclapply(df, my_mr,mc.cores = 128)
    return(processed_df)
  })
  return(processed_list)
})

combined_mr_smd_outcome <- do.call(rbind, lapply(mr_smd_outcome, function(sublist) {
  do.call(rbind, lapply(sublist, function(df1) {
    do.call(rbind, lapply(df1, function(df2) {
      return(df2)
      }))# Assuming the data frames are stored in the first element
  }))
}))
  
combined_mr_metabolite_outcome <- do.call(rbind, lapply(mr_metabolite_outcome, function(sublist) {
  do.call(rbind, lapply(sublist, function(df1) {
    do.call(rbind, lapply(df1, function(df2) {
      return(df2)
    }))# Assuming the data frames are stored in the first element
  }))
}))

combined_mr_metabolite_outcome$pval_adjusted=p.adjust(combined_mr_metabolite_outcome$pval,method = 'fdr',
                                                   n=nrow(combined_mr_metabolite_outcome))
combined_mr_smd_outcome$pval_adjusted=p.adjust(combined_mr_smd_outcome$pval,method = 'fdr',
                                                      n=nrow(combined_mr_smd_outcome))

combined_mr_metabolite_outcome2=combined_mr_metabolite_outcome[,3:10]
combined_mr_smd_outcome2=combined_mr_smd_outcome[,3:10]


get_presso=function(H_data1){
  H_data_rm_palindromic1 <- H_data1[H_data1$palindromic == F,]
  presso1=mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                    SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = H_data_rm_palindromic1, 
                    NbDistribution = 1000,  SignifThreshold = 0.05)
  sbg1=presso1$`MR-PRESSO results`$`Global Test`$RSSobs
  sbg2=presso1$`MR-PRESSO results`$`Global Test`$Pvalue
  prdf1=presso1$`Main MR results`
  prdf1$Global_test_RSSobs=sbg1
  prdf1$Global_test_Pvalue=sbg2
  prdf1$Exposure=H_data1$exposure[1:2]
  prdf1$Outcome=H_data1$outcome[1:2]
  prdf2=prdf1[,c(1,9,2:8)]
  return(prdf2)
}

presso_smd_outcome <- lapply(mental_outcome, function(sublist) {
  processed_list <- lapply(sublist, function(df) {
    processed_df <- mclapply(df, get_presso,mc.cores = 128)
    return(processed_df)
  })
  return(processed_list)
})

presso_metabolite_outcome <- lapply(metabolite_outcome, function(sublist) {
  processed_list <- lapply(sublist, function(df) {
    processed_df <- mclapply(df, get_presso,mc.cores = 128)
    return(processed_df)
  })
  return(processed_list)
})

combined_presso_metabolite_outcome <- do.call(rbind, lapply(presso_metabolite_outcome, function(sublist) {
  do.call(rbind, lapply(sublist, function(df1) {
    do.call(rbind, lapply(df1, function(df2) {
      return(df2)
    }))# Assuming the data frames are stored in the first element
  }))
}))

combined_presso_smd_outcome <- do.call(rbind, lapply(presso_smd_outcome, function(sublist) {
  do.call(rbind, lapply(sublist, function(df1) {
    do.call(rbind, lapply(df1, function(df2) {
      return(df2)
    }))# Assuming the data frames are stored in the first element
  }))
}))

combined_presso_metabolite_outcome2=combined_presso_metabolite_outcome[complete.cases(combined_presso_metabolite_outcome),]
combined_presso_smd_outcome2=combined_presso_smd_outcome[complete.cases(combined_presso_smd_outcome),]
combined_presso_metabolite_outcome2$Pval_adjusted=p.adjust(combined_presso_metabolite_outcome2$`P-value`,method = 'fdr',
                                                           n = nrow(combined_presso_metabolite_outcome2))
combined_presso_smd_outcome2$Pval_adjusted=p.adjust(combined_presso_smd_outcome2$`P-value`,method = 'fdr',
                                                           n = nrow(combined_presso_smd_outcome2))

combined_presso_metabolite_outcome3=combined_presso_metabolite_outcome2[,c(1:7,10,8,9)]
combined_presso_metabolite_outcome4=combined_presso_metabolite_outcome2[,c(1:7,10)]
combined_presso_smd_outcome3=combined_presso_smd_outcome2[,c(1:7,10,8,9)]
combined_presso_smd_outcome4=combined_presso_smd_outcome2[,c(1:7,10)]

out1=args[3]

write.table(combined_mr_metabolite_outcome2,file = paste0(out1,'_exposure_smd_outcome_metabolite.csv'),
            sep='\t',quote = F,row.names = F)
write.table(combined_mr_smd_outcome2,file = paste0(out1,'_exposure_metabolite_outcome_smd.csv'),
            sep='\t',quote = F,row.names = F)
write.table(combined_presso_metabolite_outcome3,file = paste0(out1,'_exposure_smd_outcome_metabolite_presso_global.csv'),
            sep='\t',quote = F,row.names = F)
write.table(combined_presso_metabolite_outcome4,file = paste0(out1,'_exposure_smd_outcome_metabolite_presso.csv'),
            sep='\t',quote = F,row.names = F)
write.table(combined_presso_smd_outcome3,file = paste0(out1,'_exposure_metabolite_outcome_smd_presso_global.csv'),
            sep='\t',quote = F,row.names = F)
write.table(combined_presso_smd_outcome4,file = paste0(out1,'_exposure_metabolite_outcome_smd_presso.csv'),
            sep='\t',quote = F,row.names = F)