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
library(stringr)

gland_eqtl <- fread("/home/zzb/project/neuroblastoma-MR/gwas/gland_eqtl.txt")
gland_eqtl <- gland_eqtl[,c("SNP", "A1", "A2", "b", "SE", "p", "Gene", "BP")]
colnames(gland_eqtl) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "gene", "position")
gland_eqtl$samplesize <- 233
gland_eqtl <- subset(gland_eqtl, pval < 5e-8)
# clump_data1 <- ld_clump(
#   tibble(rsid=gland_eqtl$SNP, pval=gland_eqtl$pval),
#   plink_bin = "/home/zzb/package/plink-1.07-x86_64/plink",
#   bfile = "/home/zzb/raw_data/clum_ref/1kg.v3/EUR"
# )
# saveRDS(clump_data1, "/home/zzb/project/neuroblastoma-MR/tem_res/gland_eqtl_clump.rds")
clump_data1 <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/Adrenal_Gland_clump.rds")
gland_eqtl <- subset(gland_eqtl, SNP %in% clump_data1$rsid)

gwas <- fread("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt")
colnames(gwas) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "samplesize")
out_data <- format_data(gwas, type="outcome")
fun1 <- function(sub_g){
  sub_all <- subset(gland_eqtl, gene == sub_g)
  if (dim(sub_all)[1] == 0){
    return(NA)
    stop()
  }
  exp_data <- format_data(sub_all,  type="exposure")
  harm_data <- harmonise_data(
    exposure_dat = exp_data,
    outcome_dat = out_data
  )
  if (dim(harm_data)[1] == 0){
    return(NA)
    stop()
  }
  harm_data$mr_keep="TRUE"
  harm_data$mr_keep=as.logical(harm_data$mr_keep)
  if (dim(harm_data)[1] > 1){
    test_MR <- mr(harm_data,method_list=c("mr_ivw")) #,method_list=c("mr_ivw")
    test_MR$type <- 2
    # heterogeneity <- mr_heterogeneity(harm_data)
  }
  if (dim(harm_data)[1] == 1){
    test_MR <- mr(harm_data,method_list=c("mr_wald_ratio")) #,method_list=c("mr_wald_ratio")
    test_MR$type <- 1
    
  }
  OR <-generate_odds_ratios(test_MR)
  tmp_data <- data.frame(
    gene = sub_g,
    or =  ifelse(dim(OR)[1] > 0, OR$or, NA),
    pvalue = ifelse(dim(OR)[1] > 0, OR$pval, NA),
    se = ifelse(dim(OR)[1] > 0, OR$se, NA),
    or_lci = ifelse(dim(OR)[1] > 0, OR$or_lci95, NA),
    or_uci = ifelse(dim(OR)[1] > 0, OR$or_uci95, NA),
    beta = ifelse(dim(OR)[1] > 0, OR$b, NA)
  )
  r2 <- 2*harm_data$eaf.exposure*(1-harm_data$eaf.exposure)*(harm_data$beta.exposure**2)
  r2 <- sum(r2)
  f_sta <- (3301 - dim(harm_data)[1] - 1)/(dim(harm_data)[1]) *((r2**2)/(1- r2**2))
  tmp_data$r2 <- r2
  tmp_data$f_sta <- f_sta
  return(tmp_data)
}
g1 <- unique(gland_eqtl$gene)
res <- mclapply(g1, fun1, mc.cores = 40)
names(res) <- g1
for (i in names(res)){
  tmp_data <- res[[i]]
  if(class(tmp_data) != "data.frame"){
    res[[i]] <- NULL
  }
}
res <- rbindlist(res)
res$FDR <- p.adjust(res$pvalue, method = "bonferroni")
saveRDS(res, "/home/zzb/project/neuroblastoma-MR/tem_res/gtex_ieu-a-816-prot-a-2003_clump.rds")

#####
#反向MR
gland_eqtl <- fread("/home/zzb/project/neuroblastoma-MR/gwas/gland_eqtl.txt")
gland_eqtl <- gland_eqtl[,c("SNP", "A1", "A2", "b", "SE", "p", "Gene", "BP")]
colnames(gland_eqtl) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "gene", "position")
gland_eqtl$samplesize <- 233
gwas <- fread("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt")
colnames(gwas) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "samplesize")
gwas <- subset(gwas, pval < 0.05)
exp_data <- format_data(gwas, type="exposure")
fun2 <- function(sub_g){
  sub_all <- subset(gland_eqtl, gene == sub_g)
  if (dim(sub_all)[1] == 0){
    return(NA)
    stop()
  }
  out_data <- format_data(sub_all,  type="outcome")
  harm_data <- harmonise_data(
    exposure_dat = exp_data,
    outcome_dat = out_data
  )
  if (dim(harm_data)[1] == 0){
    return(NA)
    stop()
  }
  harm_data$mr_keep="TRUE"
  harm_data$mr_keep=as.logical(harm_data$mr_keep)
  test_MR <- mr(harm_data) 
  OR <-generate_odds_ratios(test_MR)
  return(OR)
}
g2 <- unique(gland_eqtl$gene)
re_res <- mclapply(g2, fun2, mc.cores = 40)
names(re_res) <- g2
for (i in names(re_res)){
  tmp_data <- re_res[[i]]
  if(class(tmp_data) != "data.frame"){
    re_res[[i]] <- NULL
  }
}
saveRDS(re_res, "/home/zzb/project/neuroblastoma-MR/tem_res/re_gtex_ieu-a-816-prot-a-2003.rds")

mr_res <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/gtex_ieu-a-816-prot-a-2003_clump.rds")
re_res <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/re_gtex_ieu-a-816-prot-a-2003.rds")
mr_res <- subset(mr_res, FDR < 5e-8)
re_res <- re_res[mr_res$gene]
mr_g <- c()
for (i in names(re_res)){
  f <- re_res[[i]]
  f <- subset(f, pval > 0.05)
  if (dim(f)[1] > 0){
    mr_g <- c(mr_g, i)
  }
}

smr_res <- fread("/home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-8/gtexv8_nb_smr.msmr")
smr_res$FDR <- p.adjust(smr_res$p_SMR_multi, method = "BH")
smr_res <- subset(smr_res, FDR < 0.05 & p_HEIDI > 0.05)
final_g <- intersect(smr_res$Gene, mr_g)

final_mr <- subset(mr_res, gene %in% final_g)
final_smr <- subset(smr_res, Gene %in% final_g)
write.csv(final_smr[1:8,], "/home/zzb/project/neuroblastoma-MR/final_res/adrenal_final_smr.csv",
          row.names = F)
######################################
#steiger 方向验证
gland_eqtl <- fread("/home/zzb/project/neuroblastoma-MR/gwas/gland_eqtl.txt")
gland_eqtl <- gland_eqtl[,c("SNP", "A1", "A2", "b", "SE", "p", "Gene", "BP")]
colnames(gland_eqtl) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "gene", "position")
clump_data1 <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/Adrenal_Gland_clump.rds")
gland_eqtl <- subset(gland_eqtl, SNP %in% clump_data1$rsid)
gland_eqtl$samplesize <- 233
gland_eqtl <- subset(gland_eqtl, pval < 5e-8)
steiger_snp <- subset(gland_eqtl, gene %in% final_g)
exp_data <-  format_data(steiger_snp, type="exposure")
gwas <- fread("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt")
colnames(gwas) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "samplesize")
gwas <- data.frame(gwas)
gwas$samplesize <- 8182
gwas <- subset(gwas, SNP %in%exp_data$SNP)
out_data <- format_data(gwas, type="outcome")
harm_data <- harmonise_data(exposure_dat = exp_data, outcome_dat = out_data)
harm_data$mr_keep="TRUE"
harm_data$mr_keep=as.logical(harm_data$mr_keep)
steiger_sl <-steiger_filtering(harm_data)
true_sl <- subset(steiger_sl, steiger_dir == "TRUE")
true_g <- unique(true_sl$gene.exposure)

library(phenoscanner)
num <- ceiling(dim(exp_data)[1]/100)
num1 <- 1
num2 <- 100
final_scanres <- list()
for (i in c(1:num)){
  if(num2 > dim(exp_data)[1]){
    num2 <- dim(exp_data)[1]
    test_snp <- exp_data$SNP[num1:num2]
    scan_res <- phenoscanner(snpquery=test_snp)
    scan_res <- scan_res$results
    final_scanres[[i]] <- scan_res 
    break
  }
  test_snp <- exp_data$SNP[num1:num2]
  scan_res <- phenoscanner(snpquery=test_snp)
  scan_res <- scan_res$results
  final_scanres[[i]] <- scan_res 
  num1 <- num1 + 100
  num2 <- num2 + 100
}
final_scanres <- rbindlist(final_scanres)
a <- subset(steiger_snp, SNP %in% final_scanres$snp)
final_scanres <- subset(final_scanres, p < 5e-8)


###############################################################################
# 共定位做辅助验证 步骤4
immu <- fread("/home/ps/soft/smr-1.3.1-linux-x86_64/731_imm_cell/immu_p0.05.txt")
immu <- immu[,c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome",
                "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "immu")]
colnames(immu) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", 
                    "pval", "samplesize", "cell_type")
immu$cell_type <- str_split_fixed(immu$cell_type, ".vcf", 2)[,1]
immu <- immu[, c("SNP", "eaf", "pval", "samplesize", "cell_type")]
blood_eqtl <- fread("/home/zzb/raw_data/eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
blood_eqtl <- blood_eqtl[, c("SNP", "Pvalue", "GeneSymbol", "NrSamples")]
colnames(blood_eqtl) <- c("SNP", "pval", "gene", "samplesize")
input <- merge(immu, blood_eqtl, by = "SNP", all=FALSE,
               suffixes=c("_immu","_eqtl"))
input$maf <- ifelse(input$eaf > 0.5, 1-input$eaf, input$eaf)
input <- subset(input, !(is.na(input$maf)))
input <- input[order(input$pval_immu),]
input <- input[!(duplicated(input$SNP)),]
input$meta <- paste0(input$gene, "_", input$cell_type)
input <- as.data.frame(input)
PPH_results <- matrix(data = NA, nrow = length(unique(input$meta)), ncol = 6)
gene <- unique(input$meta)
for (i in 1:length(gene)){
  input2 <- input[which(input$meta == gene[i]),]
  result <- coloc.abf(dataset1=list(pvalues=input2$pval_eqtl, 
                                    type="quant", N=min(input2$samplesize_eqtl), snp = input2$SNP), 
                      dataset2=list(pvalues=input2$pval_immu, type="quant", 
                                    N=input2$samplesize_immu, snp = input2$SNP),
                      MAF = input2$maf 
  )
  PPH_results[i,1:5] <- as.numeric(result$summary[2:6])
  PPH_results[i,6] <- as.character(gene[i])
}
PPH_results <- as.data.frame(PPH_results)
colnames(PPH_results) <- c("PPH0","PPH1","PPH2","PPH3","PPH4","gene")
PPH_results$PPH4 <- as.numeric(PPH_results$PPH4)
PPH_results$PPH3 <- as.numeric(PPH_results$PPH3)
PPH_results$immu <- str_split_fixed(PPH_results$gene, "-a-", 2)[,2]
saveRDS(PPH_results, "/home/zzb/project/neuroblastoma-MR/tem_res/blood_eqtl-immu-coloc.rds")
PPH_results <- subset(PPH_results, PPH4 > 0.85)
res_coloc <- subset(input, gene %in% PPH_results$gene)
################################################################################
## 做nb与immu之间MR 步骤1
gwas <- fread("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt")
colnames(gwas) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "samplesize")
out_data <- format_data(gwas, type="outcome")
fun1 <- function(x){
  sub_all <- immu[[x]]
  sub_all$Phenotype <- x
  sub_all <- subset(sub_all, pval.outcome < 1e-5)
  sub_all$maf <- ifelse(sub_all$eaf.outcome > 0.5,
                        1-sub_all$eaf.outcome, sub_all$eaf.outcome)
  sub_all <- sub_all %>% dplyr::mutate(R2=2*(1-maf)*maf*beta.outcome**2,
                                       F.value=(R2/(1-R2))*((samplesize.outcome-1-1)/1)) %>%
    dplyr::filter(F.value>10,maf>0.01)
  if (dim(sub_all)[1] == 0){
    return(NA)
    stop()
  }
  clump_data <- readRDS(paste0("/home/zzb/project/neuroblastoma-MR/immu_cell/all_clump/",
                               x, ".rds"))
  sub_all <- subset(sub_all, SNP %in% clump_data$rsid)
  if (dim(sub_all)[1] == 0){
    return(NA)
    stop()
  }
  sub_all <- sub_all[, c("SNP", "other_allele.outcome", "effect_allele.outcome",
                         "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome",
                         "Phenotype", "samplesize.outcome")]
  colnames(sub_all) <- c("SNP", "other_allele", "effect_allele", "beta", "se",
                         "pval", "eaf", "Phenotype", "samplesize")
  exp_data <- format_data(sub_all,  type="exposure")
  harm_data <- harmonise_data(
    exposure_dat = exp_data,
    outcome_dat = out_data
  )
  if (dim(harm_data)[1] == 0){
    return(NA)
    stop()
  }
  harm_data$mr_keep="TRUE"
  harm_data$mr_keep=as.logical(harm_data$mr_keep)
  if (dim(harm_data)[1] > 1){
    test_MR <- mr(harm_data,method_list=c("mr_ivw")) #,method_list=c("mr_ivw")
    # heterogeneity <- mr_heterogeneity(harm_data)
  }
  if (dim(harm_data)[1] == 1){
    test_MR <- mr(harm_data,method_list=c("mr_wald_ratio")) #,method_list=c("mr_wald_ratio")
    
  }
  OR <-generate_odds_ratios(test_MR)
  OR$type <- x
  return(OR)
}
immu <- readRDS("/home/zzb/project/neuroblastoma-MR/immu_cell/all.rds")
g1 <- names(immu)
res1 <- mclapply(g1, fun1, mc.cores = 40)
names(res1) <- g1
saveRDS(res1, "/home/zzb/project/neuroblastoma-MR/tem_res/immu_prot-a-2003_clump_mr.rds")
#####上步反向MR 步骤5
gwas <- fread("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt")
colnames(gwas) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "samplesize")
gwas <- subset(gwas, pval < 0.05)
gwas <- data.frame(gwas)
exp_data <- format_data(gwas, type="exposure")
fun5 <- function(x){
  sub_all <- immu[[x]]
  out_data <- sub_all
  out_data$id.outcome <- x
  out_data$outcome <- "immu"
  harm_data <- harmonise_data(
    exposure_dat = exp_data,
    outcome_dat = out_data
  )
  if (dim(harm_data)[1] == 0){
    return(NA)
    stop()
  }
  harm_data$mr_keep="TRUE"
  harm_data$mr_keep=as.logical(harm_data$mr_keep)
  test_MR <- mr(harm_data) #,method_list=c("mr_wald_ratio")
  test_MR$type <- x
  return(test_MR)
}
immu <- readRDS("/home/zzb/project/neuroblastoma-MR/immu_cell/all.rds")
g5 <- names(immu)
res5 <- mclapply(g5, fun5, mc.cores = 20)
names(res5) <- g5
saveRDS(res5, "/home/zzb/project/neuroblastoma-MR/tem_res/rev_nb_immu.rds")
############################################################################
mr_res1 <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/immu_prot-a-2003_clump_mr.rds")
for (i in names(mr_res1)){
  tmp <- mr_res1[[i]]
  if(class(tmp) != "data.frame"){
    mr_res1[[i]] <- NULL
  }
}
mr_res1 <- rbindlist(mr_res1)
mr_res1$FDR <- p.adjust(mr_res1$pval, method = "bonferroni")
mr_res1 <- subset(mr_res1, FDR < 5e-8)
rev_res5 <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/rev_nb_immu.rds")
rev_res5i <- c()
for (i in names(rev_res5)){
  f <- rev_res5[[i]]
  if(!("Inverse variance weighted" %in% f$method)){
    next
  }
  f <- subset(f, method == "Inverse variance weighted")
  if(f$pval > 0.05){
    rev_res5i <- c(rev_res5i, i)
  }
}
mr_15_i <- intersect(rev_res5i, mr_res1$type)

write.table(mr_15_i, "/home/zzb/project/neuroblastoma-MR/tem_res/immu_np-mrresult.txt",
            row.names = F, quote = F, col.names = F)
#####################################################
##做血液eqtl和nb之间mr  以排除直接影响的基因 步骤3
gwas <- fread("/home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt")
colnames(gwas) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "samplesize")
out_data <- format_data(gwas, type="outcome")
blood_eqtl <- fread("/home/zzb/project/IPF-SMR/data/extract_qtl/eqtl-full-summary.txt")
ch6 <- subset(blood_eqtl, blood_eqtl$Chr == 6)
ch6 <- subset(ch6, BP >= 25119106 & BP <= 33854733)
ch8 <- subset(blood_eqtl, blood_eqtl$Chr == 8)
ch8 <- subset(ch8, BP >= 7200000 & BP <= 12500000)
mhc_snp <- c(ch6$SNP, ch8$SNP)
blood_eqtl <- subset(blood_eqtl, !(blood_eqtl$SNP %in% mhc_snp))
blood_eqtl <- subset(blood_eqtl, p < 5e-8)
blood_eqtl$ID <- "boold-eqtl"
blood_eqtl <- subset(blood_eqtl, !(is.na(blood_eqtl$Gene)))
blood_eqtl <- blood_eqtl[,c("SNP", "A1", "A2", "Freq", "b", "SE", "p", "Gene")]
colnames(blood_eqtl) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "gene")
clump_data3 <- readRDS("/home/zzb/raw_data/eqtl/smr/boold_eqtl_clump.rds")
blood_eqtl <- subset(blood_eqtl, SNP %in%clump_data3$rsid)
blood_eqtl <- data.frame(blood_eqtl)
fun3 <- function(x){
  sub_eqtl <- subset(blood_eqtl, gene == x)
  exp_data <- format_data(sub_eqtl, type="exposure")
  harm_data <- harmonise_data(
    exposure_dat = exp_data,
    outcome_dat = out_data
  )
  if (dim(harm_data)[1] == 0){
    return(NA)
    stop()
  }
  harm_data$mr_keep="TRUE"
  harm_data$mr_keep=as.logical(harm_data$mr_keep)
  test_MR <- mr(harm_data) 
  OR <-generate_odds_ratios(test_MR)
  OR$gene <- x
  return(OR)
}
g3 <- unique(blood_eqtl$gene)
res3 <- mclapply(g3, fun3, mc.cores = 60)
names(res3) <- g3
saveRDS(res3, "/home/zzb/project/neuroblastoma-MR/tem_res/eqtl_nb_clump-mr.rds")
mr_res3 <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/eqtl_nb_clump-mr.rds")
rev_res3 <- c()
for (i in names(mr_res3)){
  tmp <- mr_res3[[i]]
  if(class(tmp) != "data.frame"){
    mr_res3[[i]] <- NULL
    next
  }
  if(dim(tmp)[1] == 5){
    tmp <- subset(tmp,method == "Inverse variance weighted")
    if(tmp$pval > 0.05){
      rev_res3 <- c(rev_res3, i)
    }
  }
  if(dim(tmp)[1] == 1){
    if (is.na(tmp$pval)){
      next()
    }
    if(tmp$pval > 0.05){
      rev_res3 <- c(rev_res3, i)
    } 
  }
}
rev_res3 <- unique(rev_res3)
final_rev_mr3 <- mr_res3[rev_res3] 
final_rev_mr3 <- rbindlist(final_rev_mr3)
##################################################################
##血液与immu之间mr 步骤2
mr_15_i <- read.table("/home/zzb/project/neuroblastoma-MR/tem_res/immu_np-mrresult.txt",
                       header =F)
immu <- readRDS("/home/zzb/project/neuroblastoma-MR/immu_cell/all.rds")
immu <- immu[mr_15_i$V1]
blood_eqtl <- fread("/home/zzb/project/IPF-SMR/data/extract_qtl/eqtl-full-summary.txt")
ch6 <- subset(blood_eqtl, blood_eqtl$Chr == 6)
ch6 <- subset(ch6, BP >= 25119106 & BP <= 33854733)
ch8 <- subset(blood_eqtl, blood_eqtl$Chr == 8)
ch8 <- subset(ch8, BP >= 7200000 & BP <= 12500000)
mhc_snp <- c(ch6$SNP, ch8$SNP)
blood_eqtl <- subset(blood_eqtl, !(blood_eqtl$SNP %in% mhc_snp))
blood_eqtl <- subset(blood_eqtl, p < 5e-8)
blood_eqtl$ID <- "boold-eqtl"
blood_eqtl$maf <- ifelse(blood_eqtl$Freq > 0.5, 1-blood_eqtl$Freq, blood_eqtl$Freq)
blood_eqtl <- subset(blood_eqtl, !(is.na(blood_eqtl$Gene)))
blood_eqtl <- blood_eqtl[,c("SNP", "A1", "A2", "Freq", "b", "SE", "p", "Gene")]
colnames(blood_eqtl) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "gene")
clump_data2 <- readRDS("/home/zzb/raw_data/eqtl/smr/boold_eqtl_clump.rds")
blood_eqtl <- subset(blood_eqtl, SNP %in%clump_data2$rsid)
blood_eqtl <- data.frame(blood_eqtl)
fun2 <- function(x){
  sub_eqtl <- subset(blood_eqtl, gene == x)
  exp_data <- format_data(sub_eqtl,  type="exposure")
  res <- list()
  for (i in names(immu)){
    out_data <- immu[[i]]
    out_data$id.outcome <- i
    out_data$outcome <- "immu"
    harm_data <- harmonise_data(
      exposure_dat = exp_data,
      outcome_dat = out_data
    )
    if (dim(harm_data)[1] == 0){
      next
    }
    harm_data$mr_keep="TRUE"
    harm_data$mr_keep=as.logical(harm_data$mr_keep)
    if (dim(harm_data)[1] > 1){
      test_MR <- mr(harm_data,method_list=c("mr_ivw")) 
    }
    if (dim(harm_data)[1] == 1){
      test_MR <- mr(harm_data,method_list=c("mr_wald_ratio"))    
    }
    OR <-generate_odds_ratios(test_MR)
    OR$cell_type <- i
    res[[i]] <- OR
  }
  # res <- rbindlist(res)
  return(res)
}
g2 <- unique(blood_eqtl$gene)
res2 <- mclapply(g2, fun2, mc.cores = 60)
names(res2) <- g2
for (i in names(res2)){
  tmp <- res2[[i]]
  for (j in names(tmp)){
    tmp2 <- tmp[[j]]
    if(class(tmp2) != "data.frame"){
      tmp2[[j]] <- NULL
      next
    }
  }
  tmp <- rbindlist(tmp)
  res2[[i]] <- tmp
}
for (i in names(res2)){
  tmp <- res2[[i]]
  if(dim(tmp)[1] == 0){
    res2[[i]] <- NULL
    next
  }
  tmp$FDR <-  p.adjust(tmp$pval, method = "bonferroni")
  res2[[i]] <- tmp
}
saveRDS(res2, "/home/zzb/project/neuroblastoma-MR/tem_res/eqtl_immu_clump-mr.rds")
###############################################################################
mr_res2 <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/eqtl_immu_clump-mr.rds")
test_g <- intersect(names(mr_res2), rev_res3)
mr_res2 <- mr_res2[test_g]
for (i in names(mr_res2)){
  tmp <- mr_res2[[i]]
  tmp <- subset(tmp, FDR < 5e-8)
  if(dim(tmp)[1] == 0){
    mr_res2[[i]] <- NULL
    next
  }
}
##用这些基因去做immu和blood之间的反向MR
write.table(names(mr_res2), "/home/zzb/project/neuroblastoma-MR/tem_res/rev_immu-blood-genn.txt",
            col.names = F, row.names = F, quote = F)
#######################################################################
###immu和blood反向MR， 步骤6
test_gene6 <- read.table("/home/zzb/project/neuroblastoma-MR/tem_res/rev_immu-blood-genn.txt",
                         header = F)
test_immu6 <- read.table("/home/zzb/project/neuroblastoma-MR/tem_res/immu_np-mrresult.txt",
                         header = F)
immu <- readRDS("/home/zzb/project/neuroblastoma-MR/immu_cell/all.rds")
immu <- immu[test_immu6$V1]
blood_eqtl <- fread("/home/zzb/project/IPF-SMR/data/extract_qtl/eqtl-full-summary.txt")
blood_eqtl <- blood_eqtl[,c("SNP", "A1", "A2", "Freq", "b", "SE", "p", "Gene")]
colnames(blood_eqtl) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "gene")
blood_eqtl <- subset(blood_eqtl, gene %in% test_gene6$V1)
blood_eqtl <- data.frame(blood_eqtl)
fun6 <- function(x){
  sub_eqtl <- subset(blood_eqtl, gene == x)
  out_data <- format_data(sub_eqtl,  type="outcome")
  res <- list()
  for (i in names(immu)){
    exp_data <- immu[[i]]
    exp_data <- exp_data[,c("SNP", "effect_allele.outcome", "other_allele.outcome",
                            "beta.outcome", "se.outcome", "pval.outcome", "eaf.outcome",
                            "samplesize.outcome")]
    colnames(exp_data) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "eaf",
                            "samplesize")
    clump_data6 <- readRDS(paste0("/home/zzb/project/neuroblastoma-MR/immu_cell/all_clump/",
                                  i, ".rds"))
    exp_data <- subset(exp_data, SNP %in% clump_data6$rsid)
    if(dim(exp_data)[1] == 0){
      next
    }
    exp_data <- format_data(exp_data,  type="exposure")
    harm_data <- harmonise_data(
      exposure_dat = exp_data,
      outcome_dat = out_data
    )
    if (dim(harm_data)[1] == 0){
      next
    }
    harm_data$mr_keep="TRUE"
    harm_data$mr_keep=as.logical(harm_data$mr_keep)
    test_MR <- mr(harm_data) 
    OR <-generate_odds_ratios(test_MR)
    OR$cell_type <- i
    res[[i]] <- OR
  }
  # res <- rbindlist(res)
  return(res)
}
g6 <- unique(blood_eqtl$gene)
res6 <- mclapply(g6, fun6, mc.cores = 40)
names(res6) <- g6
saveRDS(res6, "/home/zzb/project/neuroblastoma-MR/tem_res/rev_blood_immu-mr6.rds")
#########################################################################
mr_res2 <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/eqtl_immu_clump-mr.rds")
rev_res6 <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/rev_blood_immu-mr6.rds")
mr_res2 <- mr_res2[names(rev_res6)]
test1 <- c()
test2 <- c()
for (i in names(mr_res2)){
  tmp <- mr_res2[[i]] 
  tmp <- subset(tmp, FDR < 5e-8)
  tmp$cell_type <- str_split_fixed(tmp$cell_type, ".vcf", 2)[,1]
  tmp$g_i <- paste0(i, "_", tmp$cell_type)
  test1 <- c(test1, tmp$g_i)
}
for (i in names(rev_res6)){
  tmp <- rev_res6[[i]]
  for (j in names(tmp)){
    tmp2 <- tmp[[j]]
    immu_t <- str_split_fixed(j, ".vcf", 2)[,1]
    immu_t <- paste0(i, "_", immu_t)
    if(dim(tmp2)[1] == 5){
      tmp2 <- subset(tmp2,method == "Inverse variance weighted")
      if(tmp2$pval > 0.05){
        test2 <- c(test2, immu_t)
      }
    }
    if(dim(tmp2)[1] == 1){
      if (is.na(tmp2$pval)){
        next()
      }
      if(tmp2$pval > 0.05){
        test2 <- c(test2, immu_t)
      } 
    }
  }
}
mr_26_gi <- intersect(test1, test2)
mr_26_gi <- as.data.frame(mr_26_gi)
colnames(mr_26_gi) <- "gene_immu"
mr_26_gi$gene <-str_split_fixed(mr_26_gi$gene_immu, "_", 2)[,1]
mr_26_gi$immu <-str_split_fixed(mr_26_gi$gene_immu, "_", 2)[,2] 
#############################################################################
final_mr2 <- list()
for (i in mr_26_gi$gene_immu){
  a <- str_split_fixed(i, "_", 2)[,1]
  tmp <- mr_res2[[a]]
  tmp$cell_type <- str_split_fixed(tmp$cell_type, ".vcf", 2)[,1]
  tmp$Gene <- a
  b <- str_split_fixed(i, "_", 2)[,2]
  tmp <- as.data.frame(tmp)
  tmp <- tmp[which(tmp$cell_type == b),]
  final_mr2[[i]] <- tmp
}
final_mr2 <- rbindlist(final_mr2)
test_a <- colnames(final_mr2)
test_a[15] <- "id"
colnames(final_mr2) <- test_a
final_mr2 <- merge(final_mr2, anno.data, by= "id")
write.csv(final_mr2, "/home/zzb/project/neuroblastoma-MR/final_res/immu_final_mr2.csv")
mr_res1 <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/immu_prot-a-2003_clump_mr.rds")
for (i in names(mr_res1)){
  tmp <- mr_res1[[i]]
  if(class(tmp) != "data.frame"){
    mr_res1[[i]] <- NULL
  }
}
mr_res1 <- rbindlist(mr_res1)
mr_res1$FDR <- p.adjust(mr_res1$pval, method = "bonferroni")
mr_res1 <- subset(mr_res1, FDR < 5e-8)
mr_res1$type <- str_split_fixed(mr_res1$type, ".vcf", 2)[,1]
final_mr1 <- subset(mr_res1, type %in% mr_26_gi$immu)
test_a <- colnames(final_mr1)
test_a[15] <- "id"
colnames(final_mr1) <- test_a
final_mr1 <- merge(final_mr1, anno.data, by= "id")
write.csv(final_mr1, "/home/zzb/project/neuroblastoma-MR/final_res/immu_final_mr1.csv",
          row.names = F)
#######################################################
blood_immu_coloc <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/blood_eqtl-immu-coloc.rds")
blood_immu_coloc <- subset(blood_immu_coloc, PPH4 > 0.5)
blood_immu_coloc$gene2 <- str_split_fixed(blood_immu_coloc$gene, "_", 2)[,1]
blood_immu_coloc$id <- str_split_fixed(blood_immu_coloc$gene, "_", 2)[,2]
blood_immu_coloc <- subset(blood_immu_coloc, gene %in% mr_26_gi$gene_immu)
intersect(mr_26_gi$gene_immu, blood_immu_coloc$gene)
blood_immu_coloc <- merge(blood_immu_coloc, anno.data, by = "id")
rev3 <- mr_res3[unique(mr_26_gi$gene)]
rev3 <- rbindlist(rev3)
write.csv(rev3, "/home/zzb/project/neuroblastoma-MR/final_res/immu_rev3.csv")
rev5 <- rev_res5[paste0(unique(mr_26_gi$immu), ".vcf.gz")]
rev5 <- rbindlist(rev5)
rev5$id <- str_split_fixed(rev5$type, ".vcf", 2)[,1]
anno.data <- openxlsx::read.xlsx("/home/ps/soft/smr-1.3.1-linux-x86_64/731_imm_cell/immune.cell.xlsx")
anno.data <- anno.data[, c("Panel", "Trait", "id")]
rev5 <- merge(rev5, anno.data, by= "id")
write.csv(rev5, "/home/zzb/project/neuroblastoma-MR/final_res/immu_rev5.csv")
rev6 <- list()
for (i in mr_26_gi$gene_immu){
  a <- str_split_fixed(i, "_", 2)[,1]
  b <- str_split_fixed(i, "_", 2)[,2]
  b <- paste0(b, ".vcf.gz")
  tmp <- rev_res6[[a]]
  tmp2 <- tmp[[b]]
  tmp2$gi <- i
  rev6[[i]] <- tmp2
}
rev6 <- rbindlist(rev6)
rev6$gene <-str_split_fixed(rev6$gi, "_", 2)[,1]
rev6$id <-str_split_fixed(rev6$gi, "_", 2)[,2]
rev6 <- merge(rev6, anno.data, by = "id")
write.csv(rev6, "/home/zzb/project/neuroblastoma-MR/final_res/immu_rev6.csv")

