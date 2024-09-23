library(coloc)
data(coloc_test_data)
attach(coloc_test_data)
head(D3)
S4=runsusie(D4)
library(LDlinkR)
gland_eqtl <- fread("/home/zzb/project/neuroblastoma-MR/gwas/gland_eqtl.txt")
gland_eqtl <- gland_eqtl[,c("SNP", "A1", "A2", "b", "SE", "p", "Gene")]
colnames(gland_eqtl) <- c("SNP", "effect_allele", "other_allele", "beta", "se", "pval", "gene")
gland_eqtl$samplesize <- 233


# doi:10.1038/s41588-020-0684-4
immu <- fread("/home/zzb/project/neuroblastoma-MR/immu_cell/immu_p0.05.txt")
immu <- immu[,c("SNP", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome",
            "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "immu")]
colnames(immu) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", 
                    "pval", "samplesize", "cell_type")
immu$cell_type <- str_split_fixed(immu$cell_type, ".vcf", 2)[,1]


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
# anno.data <- openxlsx::read.xlsx("/home/ps/soft/smr-1.3.1-linux-x86_64/731_imm_cell/immune.cell.xlsx")
# anno.data <- subset(anno.data, anno.data$id %in% res$protein)
# res <- merge(res, anno.data, by = "id")
res$fdr <- p.adjust(res$pvalue, method = "bonferroni")
saveRDS(res, "/home/zzb/project/neuroblastoma-MR/tem_res/gtex_ieu-a-816-prot-a-2003.rds")

plasma <- fread("/home/zzb/project/epilepsy-MR/CSF_pqtl/plasma_pqtl_uk_st10.csv")
plasma$samplesize <- 54219
plasma2 <- subset(plasma, pval > 0)
input <- merge(gwas, plasma2, by = "SNP",
               all = FALSE,
               suffixes=c("_gwas","_plasma"))
input$maf <- ifelse(input$eaf_gwas < 0.5,
                    input$eaf_gwas, 1 - input$eaf_gwas)
input <- input[order(input$pval_gwas),]
input <- input[!(duplicated(input$SNP)),]
final_g <- unique(input$gene)
PPH_results <- matrix(data = NA, nrow = length(final_g), ncol = 6)
for (i in c(1:length(final_g))){
  g <- final_g[i]
  input2 <- subset(input, gene == g)
  result <- coloc.abf(dataset1=list(pvalues=input2$pval_gwas, 
                                    type="quant", N=unique(input2$samplesize_gwas), snp = input2$SNP
  ), 
  dataset2=list(pvalues=input2$pval_plasma, type="quant", 
                N=unique(input2$samplesize_plasma),snp = input2$SNP
  ),
  MAF = input2$maf, 
  )
  PPH_results[i,1:5] <- as.numeric(result$summary[2:6])
  PPH_results[i,6] <- as.character(final_g[i])
}
PPH_results <- as.data.frame(PPH_results)
colnames(PPH_results) <- c("PPH0","PPH1","PPH2","PPH3","PPH4","gene")
PPH_results$PPH4 <- as.numeric(PPH_results$PPH4)
PPH_results$PPH3 <- as.numeric(PPH_results$PPH3)

plasma <- fread("/home/zzb/project/epilepsy-MR/CSF_pqtl/plasma_pqtl_uk_st10.csv")
plasma$samplesize <- 54219
plasma2 <- subset(plasma, pval > 0)
plasma2$chr <- str_split_fixed(plasma2$ID, ":", 6)[,1]
input <- merge(gwas, plasma2, by = "SNP",
               all = FALSE,
               suffixes=c("_gwas","_plasma"))
final_chr <- unique(input$chr)
for (i in c(1:length(final_chr))){
  g <- final_chr[i]
  input2 <- subset(input, chr == g)
  LDinfo <- LDmatrix(snps = input2$SNP, 
                     pop = "EUR", r2d = "r2", 
                     token = '5e83a011fe92', 
                     file =FALSE)
  row.names(LDinfo) <- LDinfo$RS_number
  LDinfo <- LDinfo[,-1]
  LDinfo <- as.matrix(LDinfo)
  input2 <- subset(input2, SNP %in% row.names(LDinfo))
  d1 <- list()
  d2 <- list()
  d1$beta <- input2$beta_gwas
  d1$varbeta <- input2$se_gwas
  d1$N <- unique(input2$samplesize_gwas)
  d1$MAF <- ifelse(input2$eaf_gwas < 0.5,
                   input2$eaf_gwas, 1 - input2$eaf_gwas)
  d1$LD <- LDinfo
  d1$snp <- input2$SNP
  d1$type <- "quant"
  d2$beta <- input2$beta_plasma
  d2$varbeta <- input2$se_plasma
  d2$N <- max(input2$samplesize_plasma)
  d2$MAF <- ifelse(input2$eaf_plasma < 0.5,
                   input2$eaf_plasma, 1 - input2$eaf_plasma)
  d2$LD <- LDinfo
  d2$snp <- input2$SNP
  d2$type <- "quant"
  s1 <- runsusie(d1)
  s2 <- runsusie(d2)
  if(requireNamespace("susieR",quietly=TRUE)) {
    susie.res=coloc.susie(d1, d2)
    print(susie.res$summary)
  }
}
load("/home/ps/soft/smr-1.3.1-linux-x86_64/731_imm_cell/01immune.cell.5e-08.RData")
plasma <- immu.cell.5e.08
plasma <- plasma[,c("SNP", "effect_allele.exposure", "other_allele.exposure", "maf",
                    "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", 
                    "exposure", "chr.exposure")]
colnames(plasma) <- c("SNP", "effect_allele", "other_allele", "eaf", "beta", "se", 
                      "pval", "samplesize", "gene", "chr")
