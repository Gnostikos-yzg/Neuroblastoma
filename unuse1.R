.libPaths("/home/ps/R/x86_64-pc-linux-gnu-library/4.3")
library(data.table)
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(foreach)
library(agricolae)
library(ggplot2)
library(ggrepel)
library(doParallel)
library(parallel)
library(coloc)
library(VariantAnnotation)
library(gwasglue)

# vcf <- readVcf("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/prot-a-2003.vcf.gz")
# gwas <- gwasvcf_to_TwoSampleMR(vcf, type = "outcome")
# out_data <- gwas
# # brain <- fread("/home/zzb/.synapseCache/324/128759324/female_male_combined_pQTL.txt")
# # brain <- brain[,c("SNP", "A1", "A2", "BETA", "SE", "MAF", "P", "GENE")]
# # colnames(brain) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "eaf", "pval", "gene")
# # brain$Phenotype <- "brain"
# # brain$samplesize <- 1277 
# # brain <- na.omit(brain)
# brain <- fread("/home/zzb/project/epilepsy-MR/CSF_pqtl/brain-pqtl.csv")
# fun1 <- function(sub_g){
#   sub_all <- subset(brain, protein == sub_g)
#   if (dim(sub_all)[1] == 0){
#     return(NA)
#     stop()
#   }
#   exp_data <- format_data(sub_all,  type="exposure",
#   )
#   harm_data <- harmonise_data(
#     exposure_dat = exp_data,
#     outcome_dat = out_data
#   )
#   if (dim(harm_data)[1] == 0){
#     return(NA)
#     stop()
#   }
#   harm_data$mr_keep="TRUE"
#   harm_data$mr_keep=as.logical(harm_data$mr_keep)
#   if (dim(harm_data)[1] > 1){
#     test_MR <- mr(harm_data,method_list=c("mr_ivw")) #,method_list=c("mr_ivw")
#     test_MR$type <- 2
#     # heterogeneity <- mr_heterogeneity(harm_data)
#   }
#   if (dim(harm_data)[1] == 1){
#     test_MR <- mr(harm_data,method_list=c("mr_wald_ratio")) #,method_list=c("mr_wald_ratio")
#     test_MR$type <- 1
#     
#   }
#   OR <-generate_odds_ratios(test_MR)
#   tmp_data <- data.frame(
#     protein = sub_g,
#     or =  ifelse(dim(OR)[1] > 0, OR$or, NA),
#     pvalue = ifelse(dim(OR)[1] > 0, OR$pval, NA),
#     se = ifelse(dim(OR)[1] > 0, OR$se, NA),
#     or_lci = ifelse(dim(OR)[1] > 0, OR$or_lci95, NA),
#     or_uci = ifelse(dim(OR)[1] > 0, OR$or_uci95, NA),
#     beta = ifelse(dim(OR)[1] > 0, OR$b, NA)
#   )
#   r2 <- 2*harm_data$eaf.exposure*(1-harm_data$eaf.exposure)*(harm_data$beta.exposure**2)
#   r2 <- sum(r2)
#   f_sta <- (3301 - dim(harm_data)[1] - 1)/(dim(harm_data)[1]) *((r2**2)/(1- r2**2))
#   tmp_data$r2 <- r2
#   tmp_data$f_sta <- f_sta
#   return(tmp_data)
# }
# g1 <- unique(brain$protein)
# res <- mclapply(g1, fun1, mc.cores = 40)
# names(res) <- g1
# for (i in names(res)){
#   tmp_data <- res[[i]]
#   if(class(tmp_data) != "data.frame"){
#     res[[i]] <- NULL
#   }
# }
# res <- rbindlist(res)
# # res2 <- subset(res, f_sta > 10)
# res$FDR <- p.adjust(res$pvalue, method = "BH")
# res2 <- subset(res, FDR < 0.05)
# 
# ################################################################################
# gwas2 <- gwas
# gwas2 <- gwas2[,c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome",
#                 "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")]
# colnames(gwas2) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "samplesize")
# gwas2$maf <- ifelse(gwas2$eaf > 0.5, 1-gwas2$eaf, gwas2$eaf)
# input <- merge(gwas2, brain, by = "SNP",
#                all = FALSE,
#                suffixes=c("_gwas","_brain"))
# input <- input[order(input$pval_gwas),]
# input <- input[!(duplicated(input$SNP)),]
# final_g <- unique(input$protein)
# final_g <- na.omit(final_g)
# PPH_results <- matrix(data = NA, nrow = length(final_g), ncol = 6)
# for (i in c(1:length(final_g))){
#   g <- final_g[i]
#   input2 <- subset(input, protein == g)
#   result <- coloc.abf(dataset1=list(pvalues=input2$pval_gwas, 
#                                     type="quant", N=unique(input2$samplesize_gwas), snp = input2$SNP
#   ), 
#   dataset2=list(pvalues=input2$pval_brain, type="quant", 
#                 N=unique(input2$samplesize_brain),snp = input2$SNP
#   ),
#   MAF = input2$maf, 
#   )
#   PPH_results[i,1:5] <- as.numeric(result$summary[2:6])
#   PPH_results[i,6] <- as.character(final_g[i])
# }
# PPH_results <- as.data.frame(PPH_results)
# colnames(PPH_results) <- c("PPH0","PPH1","PPH2","PPH3","PPH4","gene")
# PPH_results$PPH4 <- as.numeric(PPH_results$PPH4)
# PPH_results$PPH3 <- as.numeric(PPH_results$PPH3)
# PPH_results <- subset(PPH_results, PPH4 > 0.8)
# res2 <- subset(res, protein %in% PPH_results$gene) 
# ################################################################################
# ##plasma数据
# plasma <- fread("/home/zzb/project/epilepsy-MR/CSF_pqtl/plasma_pqtl_uk_st10.csv")
# plasma$samplesize <- 54219
# out_data <- gwas
# fun1 <- function(sub_g){
#   sub_all <- subset(plasma, gene == sub_g)
#   if (dim(sub_all)[1] == 0){
#     return(NA)
#     stop()
#   }
#   exp_data <- format_data(sub_all,  type="exposure",
#   )
#   harm_data <- harmonise_data(
#     exposure_dat = exp_data,
#     outcome_dat = out_data
#   )
#   if (dim(harm_data)[1] == 0){
#     return(NA)
#     stop()
#   }
#   harm_data$mr_keep="TRUE"
#   harm_data$mr_keep=as.logical(harm_data$mr_keep)
#   if (dim(harm_data)[1] > 1){
#     test_MR <- mr(harm_data,method_list=c("mr_ivw")) #,method_list=c("mr_ivw")
#     test_MR$type <- 2
#     # heterogeneity <- mr_heterogeneity(harm_data)
#   }
#   if (dim(harm_data)[1] == 1){
#     test_MR <- mr(harm_data,method_list=c("mr_wald_ratio")) #,method_list=c("mr_wald_ratio")
#     test_MR$type <- 1
#   }
#   OR <-generate_odds_ratios(test_MR)
#   tmp_data <- data.frame(
#     protein = sub_g,
#     or =  ifelse(dim(OR)[1] > 0, OR$or, NA),
#     pvalue = ifelse(dim(OR)[1] > 0, OR$pval, NA),
#     se = ifelse(dim(OR)[1] > 0, OR$se, NA),
#     or_lci = ifelse(dim(OR)[1] > 0, OR$or_lci95, NA),
#     or_uci = ifelse(dim(OR)[1] > 0, OR$or_uci95, NA),
#     beta = ifelse(dim(OR)[1] > 0, OR$b, NA)
#   )
#   r2 <- 2*harm_data$eaf.exposure*(1-harm_data$eaf.exposure)*(harm_data$beta.exposure**2)
#   r2 <- max(r2)
#   f_sta <- (44889 - dim(harm_data)[1] - 1)/(dim(harm_data)[1]) *((r2**2)/(1- r2**2))
#   tmp_data$r2 <- r2
#   tmp_data$f_sta <- f_sta
#   return(tmp_data)
# }
# g1 <- unique(plasma$gene)
# res <- mclapply(g1, fun1, mc.cores = 40)
# names(res) <- g1
# for (i in names(res)){
#   tmp_data <- res[[i]]
#   if(class(tmp_data) != "data.frame"){
#     res[[i]] <- NULL
#   }
# }
# res <- rbindlist(res)
# res$FDR <-p.adjust(res$pvalue, method = "BH")
# res2 <- subset(res, FDR < 0.05)
# 
# ###############################################################################
# vcf <- readVcf("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/ieu-a-816.vcf.gz")
# gwas <- gwasvcf_to_TwoSampleMR(vcf, type = "outcome")
# gwas <- gwas[,c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome",
#                 "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")]
# colnames(gwas) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
# gwas <-subset(gwas, SNP != "")
# a <- c()
# for (i in unique(gwas$SNP)){
#   if (str_detect(i, ",")){
#     a <- c(a, i)
#   }
# }
# write.table(gwas,"/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_smr.txt",
#             row.names = F, quote = F)
# sample <-  c("ebi-a-GCST004883", "ebi-a-GCST004884", "ebi-a-GCST004885",
#              "ieu-a-816", "prot-a-2003")
# res <- c()
# for (i in sample){
#   smr_res <- fread(paste0("/home/zzb/project/neuroblastoma-MR/smr_res/",
#                    i, "_Adrenal_Gland.msmr"))
#   smr_res$FDR <- p.adjust(smr_res$p_SMR, method = "BH")
#   smr_res <- subset(smr_res, p_SMR_multi < 0.05 & p_HEIDI > 0.05)
#   print(paste0(i, ": ", dim(smr_res)[1]))
# }
# 
# f <- fread("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/all.txt")
# f <- f[, c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "Freq1", "P-value")]
# colnames(f) <- c("SNP", "A1", "A2", "b", "se", "freq", "p")
# f$A1 <- toupper(f$A1)
# f$A2 <- toupper(f$A2)
# f$n <- NA  
# write.table(f, "/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/all.txt",
#             row.names = F, quote = F)  
# 
# 
# 
# library(tidyverse)
# library(data.table)
# library(TwoSampleMR)
# library(rlist)
# anno.data <- openxlsx::read.xlsx("/home/ps/soft/smr-1.3.1-linux-x86_64/731_imm_cell/immune.cell.xlsx")
# ## 按照5e-08筛选
# immu.cell <- list()
# for (x in 1:nrow(anno.data)) {
#   repeat{
#     try({immu.cell <- extract_instruments(anno.data$id[x],p1=1,p2 = 1)})
#     if(exists("immu.cell")){break}
#     Sys.sleep(0.5)
#   }
#   immu.cell <- list.append(immu.cell,immu.cell)
#   rm(immu.cell)
#   print(str_c(x," has been processed !"))
# }
# immu.cell <- do.call(rbind,immu.cell) %>% 
#   dplyr::mutate(maf=ifelse(eaf.exposure<0.5,eaf.exposure,(1-eaf.exposure)),
#                 R2=2*(1-maf)*maf*beta.exposure*beta.exposure,
#                 F.value=(R2/(1-R2))*((samplesize.exposure-1-1)/1)) 
#   # %>% 
#   # dplyr::filter(F.value>10,maf>0.01)
# saveRDS(immu.cell, "/home/zzb/project/neuroblastoma-MR/immu_cell/immu_cell.rds")
# all <- list()
# for (i in c(1:8)){
#   f <- readRDS(paste0("/home/zzb/project/neuroblastoma-MR/immu_cell/split_immu/immu_",
#                i, ".rds"))
#   all <- c(all, f)
# }


# gwas <- fread("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt")
# colnames(gwas) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "samplesize")
# out_data <- format_data(gwas, type="outcome")
# fun1 <- function(x){
#   sub_all <- immu[[x]]
#   sub_all$Phenotype <- x
#   sub_all <- subset(sub_all, pval.outcome < 1e-5)
#   sub_all$maf <- ifelse(sub_all$eaf.outcome > 0.5, 
#                         1-sub_all$eaf.outcome, sub_all$eaf.outcome)
#   sub_all <- sub_all %>% dplyr::mutate(R2=2*(1-maf)*maf*beta.outcome**2,
#                            F.value=(R2/(1-R2))*((samplesize.outcome-1-1)/1)) %>%
#               dplyr::filter(F.value>10,maf>0.01)
#   if (dim(sub_all)[1] == 0){
#     return(NA)
#     stop()
#   }
#   clump_data <- ld_clump(
#     tibble(rsid=sub_all$SNP, pval=sub_all$pval.outcome, id=sub_all$Phenotype),
#     plink_bin = "/home/zzb/package/plink-1.07-x86_64/plink",
#     bfile = "/home/zzb/raw_data/clum_ref/1kg.v3/EUR",
#     clump_r2 = 0.01
#   )
#   saveRDS(clump_data, paste0("/home/zzb/project/neuroblastoma-MR/immu_cell/clump/", x, ".rds"))
#   sub_all <- subset(sub_all, SNP %in% clump_data$rsid)
#   if (dim(sub_all)[1] == 0){
#     return(NA)
#     stop()
#   }
#   sub_all <- sub_all[, c("SNP", "other_allele.outcome", "effect_allele.outcome",
#                          "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome",
#                          "Phenotype", "samplesize.outcome")]
#   colnames(sub_all) <- c("SNP", "other_allele", "effect_allele", "beta", "se",
#                          "pval", "eaf", "Phenotype", "samplesize")
#   exp_data <- format_data(sub_all,  type="exposure")
#   harm_data <- harmonise_data(
#     exposure_dat = exp_data,
#     outcome_dat = out_data
#   )
#   if (dim(harm_data)[1] == 0){
#     return(NA)
#     stop()
#   }
#   harm_data$mr_keep="TRUE"
#   harm_data$mr_keep=as.logical(harm_data$mr_keep)
#   if (dim(harm_data)[1] > 1){
#     test_MR <- mr(harm_data,method_list=c("mr_ivw")) #,method_list=c("mr_ivw")
#     # heterogeneity <- mr_heterogeneity(harm_data)
#   }
#   if (dim(harm_data)[1] == 1){
#     test_MR <- mr(harm_data,method_list=c("mr_wald_ratio")) #,method_list=c("mr_wald_ratio")
#     
#   }
#   OR <-generate_odds_ratios(test_MR)
#   OR$type <- x
#   print(paste0(x, ": isdone!"))
#   return(OR)
# }
# immu <- readRDS("/home/zzb/project/neuroblastoma-MR/immu_cell/all.rds")
# g1 <- names(immu)
# res <- mclapply(g1, fun1, mc.cores = 40)
# saveRDS("/home/zzb/project/neuroblastoma-MR/tem_res/immu_prot-a-2003_clump.rds")

gland_eqtl <- fread("/home/zzb/project/neuroblastoma-MR/gwas/gland_eqtl.txt")
gland_eqtl <- gland_eqtl[,c("SNP", "A1", "A2", "b", "SE", "p", "Gene", "BP")]
colnames(gland_eqtl) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "gene", "position")
gland_eqtl$samplesize <- 233
gland_eqtl <- subset(gland_eqtl, pval < 5e-8)
clump_data1 <- ld_clump(
  tibble(rsid=gland_eqtl$SNP, pval=gland_eqtl$pval),
  plink_bin = "/home/zzb/package/plink-1.07-x86_64/plink",
  bfile = "/home/zzb/raw_data/clum_ref/1kg.v3/EUR"
)
saveRDS(clump_data1, "/home/zzb/project/neuroblastoma-MR/tem_res/gland_eqtl_clump.rds")  
