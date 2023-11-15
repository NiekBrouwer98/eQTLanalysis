mostsig_c1 <- filtered_ordered_length_snps_c1[1:100,]
mostsig_c2 <- filtered_ordered_length_snps_c2[1:100,]

different_eqtls <- setdiff(filtered_ordered_length_snps_c1[,1], filtered_ordered_length_snps_c2[,1])
same_eqtls <- intersect(filtered_ordered_length_snps_c1[,1], filtered_ordered_length_snps_c2[,1])
length(different_eqtls)
length(same_eqtls)
same_eqtls[1]
head(filtered_ordered_length_snps_c2)

number <- 400
genes_most_sig_c1.1 <- SNP_list_ct1[[as.character(filtered_ordered_length_snps_c2[,1][[1]])]][["gene"]]
genes_most_sig_c2.1 <- SNP_list_ct2[[as.character(filtered_ordered_length_snps_c2[,1][[1]])]][["gene"]]

corresponding_genes <- genes_most_sig_c1.1[genes_most_sig_c1.1 %in% genes_most_sig_c2.1]

snp <- as.character(filtered_ordered_length_snps_c2[,1][[1]])
corresponding_genes

length(corresponding_genes)
genes_most_sig_c2.1
gene <- corresponding_genes[[10]]

# path <- 'C:/Users/niekb/Documents/BEP datasets/simulated_data'
# genotypes <- read.delim(paste(path, '/genotypes.txt', sep = ""))
# genotypes50 <- read.delim(paste(path, '/genotypes50.txt', sep = ""))
# counts_c1 <- read.delim(paste(path, '/matrix_celltype1.txt', sep=""))
# counts_c2 <- read.delim(paste(path, '/matrix_celltype2.txt', sep=""))
# counts_c1_batches <- read.delim(paste(path, '/matrix_celltype1_batches.txt', sep=""))
# counts_c2_batches <- read.delim(paste(path, '/matrix_celltype2_batches.txt', sep=""))

snp_number_list <- genotypes[[as.character(snp)]]
snp_number_list <- cbind(rownames(genotypes), snp_number_list)
counts_list_c1 <- counts_c1[[as.character(gene)]]
counts_list_c2 <- counts_c2[[as.character(gene)]]

counts_list_c1 <- cbind(counts_c1['FID'], counts_list_c1)
counts_list_c2 <- cbind(counts_c2['FID'], counts_list_c2)

#Means cell type A
mean_list_c1 <- data.frame(matrix(nrow=10))
mean_list_c1[,1] <- snp_number_list[,1]
mean_list_c1[,2] <- c(0,0,0,0,0,0,0,0,0,0)
mean_list_c1[,3] <- c(0,0,0,0,0,0,0,0,0,0)

for (i in 1:nrow(counts_list_c1['FID'])){
  for (j in 1:nrow(mean_list_c1)){
    check <- as.character(counts_list_c1[i,1]) == mean_list_c1[j,1]
    if (check){
      mean_list_c1[j,2] = mean_list_c1[j,2] + as.numeric(counts_list_c1[i,2])
      mean_list_c1[j,3] = mean_list_c1[j,3] + 1
    }
  }
}
for (i in 1:nrow(mean_list_c1)){
  mean_list_c1[i,4] <- mean_list_c1[i,2]/mean_list_c1[i,3]
}

mean_list_c1[,5] <- snp_number_list[,2]

#Means cell type B
mean_list_c2 <- data.frame(matrix(nrow=10))
mean_list_c2[,1] <- snp_number_list[,1]
mean_list_c2[,2] <- c(0,0,0,0,0,0,0,0,0,0)
mean_list_c2[,3] <- c(0,0,0,0,0,0,0,0,0,0)

for (i in 1:nrow(counts_list_c2['FID'])){
  for (j in 1:nrow(mean_list_c2)){
    check <- as.character(counts_list_c2[i,1]) == mean_list_c2[j,1]
    if (check){
      mean_list_c2[j,2] = mean_list_c2[j,2] + as.numeric(counts_list_c2[i,2])
      mean_list_c2[j,3] = mean_list_c2[j,3] + 1
    }
  }
}
head(mean_list_c2)
for (i in 1:nrow(mean_list_c2)){
  mean_list_c2[i,4] <- mean_list_c2[i,2]/mean_list_c2[i,3]
}

mean_list_c2[,5] <- snp_number_list[,2]

# list_c1 <- list()
# for (i in 1:nrow(counts_list_c1['FID'])){
#   ID <- counts_list_c1[i,1]
#   for (j in 1:50){
#     snp <- snp_number_list[j,2]
#     ID_snp <- snp_number_list[j,1]
#     check <- (ID == ID_snp)
#     if (check){
#       list_c1 <- append(list_c1, as.character(snp))
#     }
#   }
# }
# 
# dim(counts_list_c1)
# length(list_c1)
# 
# counts_list_c1$SNP <- list_c1
# counts_list_c1 <- as.data.frame(lapply(counts_list_c1, unlist))
# 
# list_c2 <- list()
# for (i in 1:nrow(counts_list_c2['FID'])){
#   ID <- counts_list_c2[i,1]
#   for (j in 1:50){
#     snp <- snp_number_list[j,2]
#     ID_snp <- snp_number_list[j,1]
#     check <- (ID == ID_snp)
#     if (check){
#       list_c2 <- append(list_c2, as.character(snp))
#     }
#   }
# }
# 
# dim(counts_list_c2)
# length(list_c2)
# 
# counts_list_c2$SNP <- list_c2
# counts_list_c2 <- as.data.frame(lapply(counts_list_c2, unlist))
# head(snp_number_list)
# head(snp_number_list[,2])

# par(mfrow=c(1,2))
# p1 <- boxplot(counts_list_c1$counts_list_c1~counts_list_c1$SNP,ylim = c(0,40),ylab = "expression", xlab= "genotypes", main = "Cell type A")
# p2 <- boxplot(counts_list_c2$counts_list_c2~counts_list_c2$SNP,ylim = c(0,40),ylab = "expression", xlab= "genotypes", main = "Cell type B")
# mostsig_c2[number,1]
# corresponding_genes[[1]]
# genes_most_sig_c2.1[[1]]

integer_snp_list_c1 <- c()
for (i in 1:nrow(mean_list_c1)){
  integer_snp_list_c1 <- append(integer_snp_list_c1, strtoi(mean_list_c1[i,5]))
}
integer_snp_list_c2 <- c()
for (i in 1:nrow(mean_list_c2)){
  integer_snp_list_c2 <- append(integer_snp_list_c2, strtoi(mean_list_c2[i,5]))
}

library(ggplot2)
corr1 <- cor.test(x=mean_list_c1[,4], y=integer_snp_list_c1, method='spearman', exact=FALSE)
corr2 <- cor.test(x=mean_list_c2[,4], y=integer_snp_list_c2, method='spearman', exact=FALSE)
corr1.1 <- cor.test(x=mean_list_c1[,4], y=integer_snp_list_c1, method='pearson', exact=FALSE)
corr2.1 <- cor.test(x=mean_list_c2[,4], y=integer_snp_list_c2, method='pearson', exact=FALSE)

corr1[["estimate"]][["rho"]]
corr2[["estimate"]][["rho"]]
corr1.1[["estimate"]][["cor"]]
corr2.1[["estimate"]][["cor"]]

library(ggplot2)
par(mfrow=c(1,2))
p1 <- boxplot(mean_list_c1[,4]~mean_list_c1[,5], ylim=c(0,700), ylab="Expression", xlab = "Genotypes", main = paste("Cell type A. r: ",round(corr1[["estimate"]][["rho"]],2), sep = ""))
stripchart(mean_list_c1[,4]~mean_list_c1[,5], vertical = TRUE, data = mean_list_c1, 
           method = "jitter", add = TRUE, pch = 20)
p2 <- boxplot(mean_list_c2[,4]~mean_list_c2[,5], ylim=c(0,700), ylab="Expression", xlab = "Genotypes", main = paste("Cell type B. r: ",round(corr2[["estimate"]][["rho"]],2),sep = ""))
stripchart(mean_list_c2[,4]~mean_list_c2[,5], vertical = TRUE, data = mean_list_c2, 
           method = "jitter", add = TRUE, pch = 20)

