#!/bin/bash
phenotype=$1
./software/gcta \
--HEreg \
--mgrm ./data/mgrm.list \
--pheno ./data/heritability/data/${phenotype}.txt \
--qcovar ./data/heritability/data/qcovar.txt \
--covar ./data/heritability/data/covar.txt \
--out ./result/heritability/result/${phenotype}