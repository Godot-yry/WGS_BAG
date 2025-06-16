#!/bin/bash

chrom=$1
pheno=$2

/home/software/regenie \
  --step 2 \
  --chr ${chrom} \
  --bed /home/ukbFiles/Caucasian_chr${chrom} \
  --phenoFile /home/data/${pheno}_pheno.txt \
  --covarFile /home/data/${pheno}_covar.txt \
  --catCovarList male,center_Vanguard,center_SC,center_DECODE,batch \
  --maxCatLevels 27 \
  --pred /home/result/step1/step1_pred.list \
  --qt \
  --bsize 1000 \
  --threads 128 \
  --write-samples \
  --print-pheno \
  --out /home/result/step2/single/step2_chr${chrom}
