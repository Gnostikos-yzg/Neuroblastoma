.libPaths("/home/ps/R/x86_64-pc-linux-gnu-library/4.3")
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
intersects <- function (...) {
  Reduce(intersect, list(...))
}
for (slp in c(2:12)){
  cis_eqtl <- fread(paste0("/home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-",
                           slp, "/eqtl_nb_smr.msmr"))
  fdr <- p.adjust(cis_eqtl$p_SMR, method = "BH")
  cis_eqtl <- data.frame(cis_eqtl, FDR = fdr)
  cis_eqtl <- subset(cis_eqtl, FDR < 0.05, p_HEIDI > 0.05)
  
  gtex <- fread(paste0("/home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-",
                       slp, "/gtexv8_nb_smr.msmr"))
  fdr <- p.adjust(gtex$p_SMR, method = "BH")
  gtex <- data.frame(gtex, FDR = fdr)
  gtex <- subset(gtex, FDR < 0.05, p_HEIDI > 0.05)
  
  mqtl <- fread(paste0("/home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-",
                       slp, "/mqtl_nb_smr.msmr"))
  mqtl <- subset(mqtl, Gene != "NA")
  fdr <- p.adjust(mqtl$p_SMR, method = "BH")
  mqtl <- data.frame(mqtl, FDR = fdr)
  mqtl <- subset(mqtl, FDR < 0.05, p_HEIDI > 0.05)
  
  num = 1
  for (i in mqtl$Gene){
    if(str_detect(i, ";")){
      i <- strsplit(i, ";")[[1]][1]
    }
    mqtl[num, "Gene"] <- i
    num <- num + 1
  }
  
  extract_c_m <- fread(paste0("/home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-",
                              slp, "/eqtl-mqtl.msmr"))
  extract_c_m <- subset(extract_c_m, Gene != "NA")
  fdr <- p.adjust(extract_c_m$p_SMR, method = "BH")
  extract_c_m <- data.frame(extract_c_m, FDR = fdr)
  extract_c_m <- subset(extract_c_m, FDR < 0.05, p_HEIDI > 0.05)
  num = 1
  for (i in extract_c_m$Gene){
    if(str_detect(i, ";")){
      i <- strsplit(i, ";")[[1]][1]
    }
    extract_c_m[num, "Gene"] <- i
    num <- num + 1
  }
  inter_gene <- intersects(cis_eqtl$Gene, mqtl$Gene, gtex$Gene, extract_c_m$Gene)
  cat(paste0(slp, ": \n", length(inter_gene), "\n", 
               length(cis_eqtl$Gene), "\n",
               length(mqtl$Gene), "\n",
               length(gtex$Gene), "\n",
               length(extract_c_m$Gene), "\n"))
}
inter_gene <- intersects(cis_eqtl$Gene, mqtl$Gene, extract_c_m$Gene)
sub_res<- subset(gtex, gtex$Gene %in% inter_gene)

mr <- readRDS("/home/zzb/project/neuroblastoma-MR/tem_res/gtex_ieu-a-816-prot-a-2003.rds")



