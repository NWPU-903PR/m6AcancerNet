options(stringsAsFactors = FALSE)
library(dplyr)
library(Hmisc)
library(DESeq2)
library(purrr)

## Input gene annotation file
gencode.v22 <- read.delim('./gencode.v22.annotation.gff3',header = FALSE)
## Only considering protein coding genes
gencode.v22 <- gencode.v22[(gencode.v22$V3 == 'gene'), ]
gene.type <- sub('^.*gene_type=(.*?);.*$', '\\1', gencode.v22$V9)
gencode.v22 <- gencode.v22[gene.type == 'protein_coding', ]
ensg.id <- sub('^.*ID=(.*?);.*$', '\\1', gencode.v22$V9)
ensg.id <- str_sub(ensg.id,1,15)
genesembal <- sub('^.*gene_name=(.*?);.*$', '\\1', gencode.v22$V9)

## Input gene expression reads count data for Tumor
exp1 <- read.table("LAML.GeneExp.HTseq.txt",header = TRUE)

## Data per-processing
exp1 <- exp1[grep('^__.*', exp1$Tags) * -1, ]
rownames(exp1) <- exp1$Tags
col.factor <- estimateSizeFactorsForMatrix(exp1[, -1])
exp <- as.data.frame(t(t(exp1[, -1]) * (1 / col.factor)))
rownames(exp) <- str_sub(rownames(exp),1,15)
exp <- exp[rownames(exp) %in% ensg.id, ]
gene.name <- genesembal[match(rownames(exp), ensg.id)]
repetition <- as.data.frame(table(gene.name))
repetition <- repetition[repetition$Freq > 1, ]
del <- c()
for (i in 1:dim(repetition)[1]) {
  rownum <- which(gene.name == repetition[i, 1])
  exp[rownum[1], ] <- colMedians(as.matrix(exp[rownum, ]))
  del <- c(del, rownum[2:length(rownum)] * -1)
}
exp <- exp[del, ]
gene.name <- gene.name[del]
rownames(exp) <- gene.name

## Removing low expressed genes
exp <- exp[rowSums(exp > 30) >= (0.8 * dim(exp)[2]), ]

saveRDS(exp,file = 'LAML_Count_Normalization.rds')


## Calculation PCC matrix
pcc <- rcorr(t(as.matrix(exp)), type = c("pearson"))
pcc.sig <- pcc$r * (pcc$P < 0.01)
diag(pcc.sig) <- 0

#convert pcc into 0/1
pcc.bin <- ifelse(near(pcc.sig,0),0,1)
