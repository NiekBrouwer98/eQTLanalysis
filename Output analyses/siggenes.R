library(DT)
summary_c1 <- as.data.frame(summary(significant_list_c1))
length_list <- list()
for (i in 1:length(significant_list_c1)){
  length_list <- append(length_list, list(length(significant_list_c1[[i]][["name"]])))
}
summary_c1 <- summary_c1[1:748,]
summary_c1[['snps']] <- length_list
summary_c1 <- cbind(summary_c1['Var1'], summary_c1['snps'])
datatable(summary_c1, colnames = c('Gene', 'siginificant SNPs'))

summary_c2 <- as.data.frame(summary(significant_list_c2))
length_list <- list()
for (i in 1:length(significant_list_c2)){
  length_list <- append(length_list, list(length(significant_list_c2[[i]][["name"]])))
}
summary_c2 <- summary_c2[1:735,]
summary_c2[['snps']] <- length_list
summary_c2 <- cbind(summary_c2['Var1'], summary_c2['snps'])
datatable(summary_c2, colnames = c('Gene', 'siginificant SNPs'))

summary_c1_batches <- as.data.frame(summary(significant_list_c1_batches))
length_list <- list()
for (i in 1:length(significant_list_c1_batches)){
  length_list <- append(length_list, list(length(significant_list_c1_batches[[i]][["name"]])))
}
summary_c1_batches <- summary_c1_batches[1:45,]
summary_c1_batches[['snps']] <- length_list
summary_c1_batches <- cbind(summary_c1_batches['Var1'], summary_c1_batches['snps'])
datatable(summary_c1_batches, colnames = c('Gene', 'siginificant SNPs'))

summary_c2_batches <- as.data.frame(summary(significant_list_c2_batches))
length_list <- list()
for (i in 1:length(significant_list_c2_batches)){
  length_list <- append(length_list, list(length(significant_list_c2_batches[[i]][["name"]])))
}
summary_c2_batches <- summary_c2_batches[1:45,]
summary_c2_batches[['snps']] <- length_list
summary_c2_batches <- cbind(summary_c2_batches['Var1'], summary_c2_batches['snps'])
datatable(summary_c2_batches, colnames = c('Gene', 'siginificant SNPs'))
