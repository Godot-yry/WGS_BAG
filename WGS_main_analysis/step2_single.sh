#!/bin/bash

chrom=$1
pheno=$2

./software/regenie \
  --step 2 \
  --chr ${chrom} \
  --bed ../Caucasian_chr${chrom} \
  --phenoFile ./data/${pheno}_pheno.txt \
  --covarFile ./data/${pheno}_covar.txt \
  --catCovarList male,center_Vanguard,center_SC,center_DECODE,batch \
  --maxCatLevels 27 \
  --pred ./result/step1/step1_pred.list \
  --qt \
  --bsize 1000 \
  --threads 128 \
  --write-samples \
  --print-pheno \
  --out ./result/step2/single/step2_chr${chrom}
