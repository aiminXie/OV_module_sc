setwd('/home/data/t060324/jobs/aimin/')
rm(list=ls())
gc()

##处理表达数据
OV_FPKM <- read.table('./data/RawData/TCGA-OV.htseq_fpkm-uq.tsv.gz',header = T,sep = '\t')

geneID <- gsub(pattern = '\\.[0-9]+$','',OV_FPKM$Ensembl_ID)
OV_FPKM <- as.matrix(OV_FPKM[,-1])
rownames(OV_FPKM) <- geneID
colnames(OV_FPKM) <- gsub(pattern = '\\.',replacement = '-',colnames(OV_FPKM))
save(OV_FPKM,file = './data/OV/OV_FPKM.RData')

OV_FPKM <- (2^OV_FPKM)-1

## FPKM2TPM
fpkm2tpm = function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

OV_TPM = apply(OV_FPKM, 2, fpkm2tpm)
colnames(OV_TPM) <- gsub(pattern = '\\.',replacement = '-',colnames(OV_TPM))
save(OV_TPM,file = './data/bulk/OV_TPM.RData')