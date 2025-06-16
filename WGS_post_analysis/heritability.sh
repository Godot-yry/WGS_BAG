#!/bin/bash
phenotype=$1
/home/software/gcta \
--HEreg \
--mgrm /home/data/mgrm.list \
--pheno /home/data/heritability/data/${phenotype}.txt \
--qcovar /home/data/heritability/data/qcovar.txt \
--covar /home/data/heritability/data/covar.txt \
--out /home/result/heritability/result/${phenotype}