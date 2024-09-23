SMR=/home/ps/soft/smr-1.3.1-linux-x86_64/smr-1.3.1
cd /home/zzb/project/neuroblastoma-MR/draw_smr/
for i in $(cat ./draw_EN.txt)
do        
    ${SMR} \
      --bfile /home/zzb/raw_data/g1000_eur/g1000_eur \
      --gwas-summary /home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt \
      --beqtl-summary /home/ps/soft/smr-1.3.1-linux-x86_64/GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite/Adrenal_Gland.lite \
      --out ./${i} \
      --plot --probe ${i} --probe-wind 500 \
      --diff-freq-prop 0.9 \
      --diff-freq 1 \
      --gene-list /home/zzb/raw_data/smr_refe_gene/glist_hg19_strand.list
done