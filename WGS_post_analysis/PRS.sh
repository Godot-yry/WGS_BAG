#!/bin/bash

pheno=$1

for chr in {1..22}
do
Rscript /home/software/PRSice/PRSice.R \
 --prsice /home/software/PRSice/PRSice_linux \
 --base /home/data/PRS/${pheno}_chr${chr}.csv \
 --snp ID \
 --chr CHROM \
 --bp GENPOS \
 --a1 ALLELE1 \
 --a2 ALLELE0 \
 --beta \
 --stat BETA \
 --pvalue P \
 --clump-kb 500kb \
 --clump-p 1.000000 \
 --clump-r2 0.100000 \
 --fastscore  \
 --bar-levels 0.00000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01 \
 --no-regress  \
 --target /home/ukbFiles/Caucasian_chr${chrom} \
 --binary-target T \
 --out /home/result/PRS/${pheno}_chr${chr}
done