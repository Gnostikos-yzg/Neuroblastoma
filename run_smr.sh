SMR=/home/ps/soft/smr-1.3.1-linux-x86_64/smr-1.3.1
for i in {2,3,4,5,6,7,8,9,10,11,12}
do
    mkdir /home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-${i}
    ${SMR} \
      --bfile /home/zzb/raw_data/g1000_eur/g1000_eur \
      --gwas-summary /home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt \
      --beqtl-summary /home/zzb/raw_data/mqtl_blood/multiple \
      --out /home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-${i}/mqtl_nb_smr \
      --diff-freq-prop 0.9 \
      --diff-freq 1 \
      --peqtl-smr 5e-0${i} \
      --smr-multi \
      --thread-num 10
    ${SMR} \
      --bfile /home/zzb/raw_data/g1000_eur/g1000_eur \
      --gwas-summary /home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt \
      --beqtl-summary /home/ps/soft/smr-1.3.1-linux-x86_64/GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite/Adrenal_Gland.lite \
      --out /home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-${i}/gtexv8_nb_smr \
      --diff-freq-prop 0.9 \
      --diff-freq 1 \
      --peqtl-smr 5e-0${i} \
      --smr-multi \
      --thread-num 10
    ${SMR} \
      --bfile /home/zzb/raw_data/g1000_eur/g1000_eur \
      --gwas-summary /home/zzb/project/neuroblastoma-MR/gwas/opengwas/smr/ieu-a-816_prot-a-20031.txt \
      --beqtl-summary /home/zzb/raw_data/eqtl/smr/cis-eQTLs-full \
      --out /home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-${i}/eqtl_nb_smr \
      --diff-freq-prop 0.9 \
      --diff-freq 1 \
      --peqtl-smr 5e-0${i} \
      --smr-multi \
      --thread-num 10
    # ${SMR} \
    #   --bfile /home/zzb/raw_data/g1000_eur/g1000_eur \
    #   --gwas-summary /home/zzb/project/IPF-SMR/data/extract_qtl/cis-eqtl-full-summary_smr.txt \
    #   --beqtl-summary  /home/zzb/raw_data/mqtl_blood/multiple \
    #   --out /home/zzb/project/neuroblastoma-MR/smr_res/smr_res_e-${i}/eqtl-mqtl \
    #   --diff-freq-prop 0.9 \
    #   --diff-freq 1 \
    #   --peqtl-smr 5e-0${i} \
    #   --smr-multi \
    #   --thread-num 10
done