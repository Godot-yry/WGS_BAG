if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(BiocManager)
library(clusterProfiler)
library(org.Hs.eg.db)
rm(list =ls())

phenos <- c("Adrenal","Artery","Brain","Esophagus","Heart","Immune","Intestine","Kidney","Liver","Lung","Muscle","Pancreas","Pituitary","Salivary","Skin","Stomach","Thyroid")
setwd('/home/data/enrichment/')

i = 1
for (i in 1:17){
  ###### READ GENE ######
  gene <- read.table(paste0('./',phenos[i],'.csv'), header=F)
  gene <- gene$V1
  
  ###### GET ENTREZID ######
  gene.df <- bitr(gene, fromType="SYMBOL",toType="ENTREZID", OrgDb = "org.Hs.eg.db")
  
  ###### GO ######
  enrich.go <- enrichGO(gene = gene.df$ENTREZID,
                        OrgDb = 'org.Hs.eg.db',
                        ont = 'ALL',
                        pAdjustMethod = 'fdr',
                        pvalueCutoff = 1,
                        qvalueCutoff = 1,
                        readable = T)
  enrich.go <- data.frame(enrich.go)
  out <- paste0('/home/result/enrichment/',phenos[i],'_enrich_GO.csv')
  write.csv(enrich.go,out,row.names = F)
  
  ###### KEGG ######
  enrich.kegg <- enrichKEGG(gene = gene.df$ENTREZID,
                            organism = 'hsa',
                            keyType = 'kegg',
                            pAdjustMethod = 'fdr',
                            pvalueCutoff = 1,
                            qvalueCutoff = 1)
  enrich.kegg <- setReadable(enrich.kegg,
                             OrgDb = "org.Hs.eg.db",
                             keyType = "ENTREZID")
  enrich.kegg <- data.frame(enrich.kegg)
  out <- paste0('/home/result/enrichment/',phenos[i],'_enrich_KEGG.csv')
  write.csv(enrich.kegg,out,row.names = F)
}



