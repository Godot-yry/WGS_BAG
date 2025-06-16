#!/bin/bash

chrom=$1
pheno=$2
maskName=$3

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
  --firth --approx \
  --anno-file /home/data/chr${chrom}/Main/${maskName}_chr${chrom}.txt \
  --set-list /home/data/chr${chrom}/Main/chr${chrom}_${maskName}.setlist \
  --mask-def /home/data/Mask/Mask_${maskName}.txt \
  --vc-tests skato,acato-full \
  --vc-maxAAF 0.01 \
  --aaf-bins 0.01,0.001,0.0001,0.00001 \ 
  --rgc-gene-p \ 
  --write-mask \
  --bsize 200 \
  --out /home/result/step2/genebased/step2_chr${chrom}_${maskName}

