#!/bin/bash
phenotype=$1
i=$2
for j in {2..6}; do
./software/gcta \
--HEreg-bivar ${i} ${j} \
--mgrm ./data/mgrm.list \
--pheno ./data/genetic_correlation/${phenotype}.txt \
--out ./result/genetic_correlation/${phenotype}_${i}_${j} \
--thread-num 95
done