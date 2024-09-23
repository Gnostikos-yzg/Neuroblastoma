library(data.table)
library(grid)
library(forestploter)
######################################################################
##绘制 肾上腺结果图
mr_res <- read.csv("/home/zzb/project/neuroblastoma-MR/final_res/adrenal_final_mr.csv",
                   header = T)
mr_res <- mr_res[, c("gene", "or", "FDR", "or_lci", "or_uci")]
colnames(mr_res) <- c("Gene", "or", "MR_P", "or_lci", "or_uci")
mr_res$'OR (95%CI)' <- paste0(round(mr_res$or, 2), "(", round(mr_res$or_lci, 2),
                              ", ", round(mr_res$or_uci, 2), ")")
# mr_res$MR_P <- round(mr_res$MR_P, 2)
smr_res <- read.csv("/home/zzb/project/neuroblastoma-MR/final_res/adrenal_final_smr.csv",
                    header = T)
write.table(smr_res$probeID, 
            "/home/zzb/project/neuroblastoma-MR/draw_smr/draw_EN.txt",
            col.names = F, row.names = F, quote = F)
smr_res <- smr_res[,c("Gene", "FDR", "p_HEIDI")]
colnames(smr_res) <- c("Gene", "SMR_P", "SMR_HEIDI")
smr_res$SMR_HEIDI <- round(smr_res$SMR_HEIDI, 2)
rev_res <- read.csv("/home/zzb/project/neuroblastoma-MR/final_res/adrenal_rev_MR.csv",
                    header = T)
rev_res <- rev_res[, c("Gene", "pval")]
rev_res$pval <- round(rev_res$pval, 2)
colnames(rev_res) <- c("Gene", "Bidirectional\nMR_P")
draw_data <- merge(mr_res, smr_res, by = "Gene")
draw_data <- merge(draw_data, rev_res, by = "Gene")
steiger <- read.csv("/home/zzb/project/neuroblastoma-MR/final_res/adrenal_steiger.csv",
                    header = T)
steiger <- steiger[, c("gene.exposure", "steiger_pval")]
colnames(steiger) <- c("Gene", "p")
steiger$'Steiger filtering' <- paste0("Passed", "(", steiger$p, ")")
steiger <- steiger[,c(1,3)]
draw_data <- merge(draw_data, steiger, by = "Gene")
draw_data$'' <- paste(rep(" ", 8), collapse = " ")
draw_data$MR_P <- as.character(draw_data$MR_P)
tm <- forest_theme(base_size = 12,  #文本的大小
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,   
                   ci_col = "red",    
                   ci_fill = "black",    
                   ci_alpha = 0.8,       
                   ci_lty = 1,            
                   ci_lwd = 1,         
                   ci_Theight = 0.1,
                   refline_lwd = 1,     
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   xlab_adjust = "center"
)
pdf("/home/zzb/project/neuroblastoma-MR/final_res/adrenal_table1.pdf",
    width=10,height=4)
p <- forest(draw_data[,c(1,11, 6, 3, 7, 8, 9, 10)],
            est = draw_data$or,
            lower = draw_data$or_lci,
            upper = draw_data$or_uci,
            sizes = 0.2,
            ci_column = 2,
            ref_line = 1,
            arrow_lab = c("Lower", "Higher"),
            xlim = c(0.80, 1.2),
            ticks_at = c(0.85, 1, 1.1),
            theme = tm
)
p <- add_border(p, part = "header")
p <- edit_plot(p, which = "background",
               gp = gpar(fill = "white"))
print(p)
dev.off()
################################################################################
###
source("/home/zzb/project/IPF-SMR/src/R/plot_SMR.r")
for (i in mr_res$Gene){
  smrdata <- ReadSMRData(paste0("/home/zzb/project/neuroblastoma-MR/draw_smr/plot/",
                               i, ".txt"))
  smrdata$SMR <- subset(smrdata$SMR, V4  == i )
  smrdata$SMR <- smrdata$SMR[order(smrdata$SMR$V8, smrdata$SMR$V9),]
  smrdata$SMR <- smrdata$SMR[!(duplicated(smrdata$SMR$V4)),]
  pdf(paste0("/home/zzb/project/neuroblastoma-MR/draw_smr/",
             i, "_locus.pdf"),width = 20, height = 10)
  p1 <- SMRLocusPlot(data=smrdata,smr_thresh=1,anno_selfdef=TRUE, my_color="steelblue2",
                     heidi_thresh = 0.05, plotWindow = 1000)
  print(p1)
  dev.off()
  pdf(paste0("/home/zzb/project/neuroblastoma-MR/draw_smr/",
            i, "_effect.pdf"), width = 15, height = 10)
  p2 <- SMREffectPlot(data=smrdata)
  print(p2)
  dev.off()
}


################################################################################
##绘制免疫细胞结果
mr2 <- read.csv("/home/zzb/project/neuroblastoma-MR/final_res/immu_final_mr2.csv",
                header = T)
mr2 <- mr2[,c("Gene", "or", "or_lci95", "or_uci95", "FDR", "cell_type")]
colnames(mr2) <- c("Gene", "or1", "or_lci1", "or_uci1", "MR1_P", "id")
anno.data <- openxlsx::read.xlsx("/home/ps/soft/smr-1.3.1-linux-x86_64/731_imm_cell/immune.cell.xlsx")
anno.data <- anno.data[, c("Panel", "Trait", "id")]
mr2 <- merge(mr2, anno.data, by = "id")
mr2$'OR1 (95%CI)' <- paste0(round(mr2$or1, 2), "(", round(mr2$or_lci1, 2),
                           ", ", round(mr2$or_uci1, 2), ")")
mr1 <- read.csv("/home/zzb/project/neuroblastoma-MR/final_res/immu_final_mr1.csv",
                header = T) 
mr1 <- mr1[,c("or", "or_lci95", "or_uci95", "type", "FDR")]
colnames(mr1) <- c("or2", "or_lci2", "or_uci2", "id", "MR2_P")
mr1$'OR2 (95%CI)' <- paste0(round(mr1$or2, 2), "(", round(mr1$or_lci2, 2),
                           ", ", round(mr1$or_uci2, 2), ")")
draw_data <- merge(mr2, mr1, by = "id")
# write.csv(draw_data, "/home/zzb/project/neuroblastoma-MR/final_res/immu_draw_OR.csv")
draw_data <- read.csv("/home/zzb/project/neuroblastoma-MR/final_res/immu_draw_OR.csv",
                      header = T, check.names = F)
draw_data$'' <- paste(rep(" ", 73), collapse = " ")
draw_data$'' <- paste(rep(" ", 73), collapse = " ")
tm <- forest_theme(base_size = 12,  #文本的大小
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,   
                   ci_col = "red",    
                   ci_fill = "black",    
                   ci_alpha = 0.8,       
                   ci_lty = 1,            
                   ci_lwd = 1,         
                   ci_Theight = 0.1,
                   refline_lwd = 1,     
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   xlab_adjust = "center"
)
for (i in unique(draw_data$Panel)){
  draw_data2 <- subset(draw_data, Panel == i)
  if(dim(draw_data2)[1] > 15){
    h <- 13
  }
  else{
    h <- 8
  }
  pdf(paste0("/home/zzb/project/neuroblastoma-MR/final_res/",
       i, "_table.pdf"),
      width=25,height=h)
  # draw_data <- draw_data[1:20,]
  p <- forest(draw_data2[,c(1,15, 5, 6, 8, 9,
                           16, 13, 14)],
              est = list(draw_data2$or1,
                         draw_data2$or2),
              lower = list(draw_data2$or_lci1,
                           draw_data2$or_lci2),
              upper = list(draw_data2$or_uci1,
                           draw_data2$or_uci2),
              sizes = 0.2,
              ci_column = c(2, 7),
              ref_line = 1,
              arrow_lab = c("Lower", "Higher"),
              xlim = list(c(0.1, 5),
                          c(0.5, 1.8)),
              ticks_at = list(c(0.5, 1, 2.5,5),
                              c(0.6, 0.8, 1, 1.2, 1.4)),
              theme = tm
  )
  p <- add_border(p, part = "header")
  p <- edit_plot(p, which = "background",
                 gp = gpar(fill = "white"))
  print(p)
  dev.off()
}
################################################################################









