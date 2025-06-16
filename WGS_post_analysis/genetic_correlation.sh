#!/bin/bash
phenotype=$1
i=$2
for j in {2..6}; do
/home/software/gcta \
--HEreg-bivar ${i} ${j} \
--mgrm /home/data/mgrm.list \
--pheno /home/data/genetic_correlation/${phenotype}.txt \
--out /home/result/genetic_correlation/${phenotype}_${i}_${j} \
--thread-num 95
done