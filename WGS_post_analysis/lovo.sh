#!/bin/bash
#SBATCH -J chr1
#SBATCH -o job%j.out
#SBATCH -e job%j.err

chr=$1
pheno=$2
maskid=$3
mask=$4
gene=$5
mafbin=$6
domain=$7

./software/regenie \
--step 2 \
--chr ${chr} \
--bed ../Caucasian_chr${chr} \
--phenoFile ./data/${pheno}_pheno.txt \
--covarFile ./data/${pheno}_covar.txt \
--phenoCol ${pheno} \
--catCovarList male,center_Vanguard,center_SC,center_DECODE,batch \
--maxCatLevels 27 \
--pred ./result/step1/step1_pred.list \
--qt \
--anno-file ./LOVO/lovofiles/${mask}_chr${chr}.txt \
--set-list ./LOVO/lovofiles/chr${chr}_${mask}.setlist \
--mask-def ./data/Mask/Mask_${maskName}.txt \
--mask-lovo ${gene}_${domain},${maskid},${mafbin} \
--aaf-bins 0.01,0.001,0.0001,0.00001 \
--bsize 1000 \
--out ./result/lovo/lovo__${pheno}__chr${chr}__${gene}__${domain}__${mask}__${mafbin}