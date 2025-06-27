#!/bin/bash

pheno=$1

./sofetware/regenie \
  --step 1 \
  --bed ../hg38 \
  --extract ./data/qc_pass.snplist \
  --phenoFile ./data/${pheno}_pheno.txt \
  --covarFile ./data/${pheno}_covar.txt \
  --catCovarList male \
  --qt \
  --bsize 1000 \
  --threads 150 \
  --lowmem \
  --lowmem-prefix  ./result/step1/tmp \
  --out ./result/step1/step1

