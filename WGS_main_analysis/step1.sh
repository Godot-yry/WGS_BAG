#!/bin/bash

pheno=$1

regenie \
  --step 1 \
  --bed /home/ukbFiles/hg38 \
  --extract /home/data/qc_pass.snplist \
  --phenoFile /home/data/${pheno}_pheno.txt \
  --covarFile /home/data/${pheno}_covar.txt \
  --catCovarList male \
  --qt \
  --bsize 1000 \
  --threads 150 \
  --lowmem \
  --lowmem-prefix  /home/result/step1/tmp \
  --out /home/result/step1/step1

