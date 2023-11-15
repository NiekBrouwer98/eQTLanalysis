if(!require(devtools)) install.packages("devtools")
devtools::install_github("XuegongLab/SCeQTL")

library(SCeQTL)
data(test)
dim(gene)
dim(snp)
checkdist(gene)
normalize(gene)
test_result <- cal.pvalue(gene, snp)
check.sample(gene[1,], snp[1,])


genotypes_matrix <- (genotypes[["genotypes"]])
colnames(genotypes_matrix)<- genotypes[["id_snps"]]
rownames(genotypes_matrix)<- genotypes[["id_samples"]]
checkdist(output_data_celltype1)
gene_matrix_celltype1 <- count_matrix_celltype1
colnames(gene_matrix_celltype1) <- names_celltype1[[1]]

dim(gene_matrix_celltype1)
dim(genotypes_matrix)

normalize(gene_matrix_celltype1)
celltype1_result <- cal.pvalue(gene_matrix_celltype1[1:10,], genotypes_matrix)
check.sample(gene_matrix_celltype1[1,], genotypes_matrix[1,])
