library(ggplot2)
library("gplots")
library(cluster)
library(plotly)
library(lattice)
library(Hmisc)
library(edgeR)
library(corrplot)
library(miRcomp)
library(miRLAB)
library(genefilter)
library(dplyr)
require(Hmisc)
require(dplyr)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(reshape2)

##Genes expression
genes_up <- read.table (%genes filter % , header = F) 
mRNA <- read.table (%rpkm genes % , header = T, sep = "\t", row.names= 1)

mRNA <- mRNA[rownames(mRNA) %in% as.vector(genes_up$V1),]

my_gene_DE <- read.table (%genes file% , header = T, sep = "\t")
my_gene_DE <- my_gene_DE[my_gene_DE$gene %in% as.vector(genes_up$V1),]

rownames(my_gene_DE) <- my_gene_DE$gene
my_gene_DE$gene <- NULL

write.table(my_gene_DE, file = %name file%, row.names=TRUE, col.names=TRUE, sep = "\t ", quote = FALSE)

##miRNA expression
miRNA <- read.table ( %rpkm /miRNA %, header = T, sep = "\t", row.names= 1)
my_miRNA_DE <- read.table (%miRNA filter Fold change and up%, header = T, sep = "\t")

##annotation file
annot <- read.delim(%annotation file%,  header = TRUE, sep = "\t")

##convert to matrix
mRNA <- as.matrix(mRNA)
miRNA <- as.matrix(miRNA)

##cmp estimate
mRNA_cpm <- (mRNA + 0.5)
miRNA_cpm <- (miRNA + 0.5)

##filter
temp <- row.names(mRNA_cpm) %in% row.names(my_gene_DE)
mRNA_cpm_DE <- subset(mRNA_cpm, temp)
mRNA_cpm_DE <- as.data.frame(mRNA_cpm_DE)

####### Subset miRNA_DE
temp <- row.names(miRNA_cpm) %in% my_miRNA_DE$Name
miRNA_cpm_DE <- subset(miRNA_cpm, temp)
miRNA_cpm_DE <- as.data.frame(miRNA_cpm_DE)

######logFC value
logFC <- function(data){
  log2_results <- row.names(data)
  for(i in 1:(ncol(data)/2)){
    tumor <- data[,paste("TUMOR_",i, sep = "")]
    normal <- data[,paste("NORMAL_",i, sep = "")]
    log2_results <- cbind(log2_results, log2(tumor/normal))
  }
  return(log2_results)
}
logFC_mRNA <- logFC(mRNA_cpm_DE)
write.table(logFC_mRNA, file = "logFC_mRNA_DE.txt", quote = F, col.names = F, row.names =F, sep= "\t")
logFC_mRNA<- read.table("logFC_mRNA_DE.txt", dec =".", sep = "\t")

logFC_miRNA <- logFC(miRNA_cpm_DE)
write.table(logFC_miRNA, file= "logFC_miRNA.txt", quote = F, col.names = F, row.names= F, sep= "\t")
logFC_miRNA<- read.table("logFC_miRNA.txt", dec =".", sep = "\t")

t_logFC_miRNA <- t(logFC_miRNA)
t_logFC_mRNA_DE <- t(logFC_mRNA)
write.table(t_logFC_miRNA, file= "t_logFC_miRNA.txt", quote = F, col.names = F, sep= "\t")
write.table(t_logFC_mRNA_DE, file = "t_logFC_mRNA_DE.txt", quote = F, col.names = F, sep= "\t")

t_logFC_mRNA_DE <- read.table("t_logFC_mRNA_DE.txt", dec =".", sep = "\t", head = TRUE, row.names = 1)
t_logFC_miRNA <- read.table("t_logFC_miRNA.txt", dec =".", sep = "\t", head = TRUE,row.names = 1)

CorrMatrix <- rcorr(as.matrix(t_logFC_miRNA),as.matrix(t_logFC_mRNA_DE), type="pearson")

Corr_r <- CorrMatrix$r
Corr_pvalue <- CorrMatrix$P

namesMiRNA <- as.vector(colnames(t_logFC_miRNA))
namesGenes <- colnames(t_logFC_mRNA_DE)

# select factor column 
tmp <- subset.matrix(Corr_pvalue, select= namesMiRNA) 
tmp2 <- subset.matrix(t(tmp), select= namesGenes)
matrix_Corr_value <- t(tmp2)

#select factor column 
tmp <- subset.matrix(Corr_r, select= namesMiRNA) 
#select genes column
tmp2 <- subset.matrix(t(tmp), select= namesGenes) 
matrix_Corr_r <- t(tmp2)

rordered_matrix_Corr_r <- melt(matrix_Corr_r)
colnames(rordered_matrix_Corr_r) <- c("Gene", "miRNA", "Correlation")
rordered_matrix_Corr_value <- melt(matrix_Corr_value)
colnames(rordered_matrix_Corr_value) <- c("Gene", "miRNA", "p_value")

final_matrix <- merge (rordered_matrix_Corr_r, rordered_matrix_Corr_value, by = c("Gene", "miRNA"))

filter_matrix <- dplyr::filter(final_matrix, final_matrix$p_value<0.05 & final_matrix$Correlation<=(-0.75))

write.table(filter_matrix[1:2], file = %name file%, row.names=FALSE, col.names=FALSE, sep = "\t", quote = FALSE)
write.table(filter_matrix, file = %name file%, row.names=FALSE, col.names=FALSE, sep = "\t", quote = FALSE)


