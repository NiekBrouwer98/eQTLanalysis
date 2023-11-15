library(GWASTools)

path = "C:/Users/niekb/plink_workspace/eqtl_misclass10p/"

SNP_list_misclass10p <- list()
for (j in 1:1500){
  assoc_table <- read.table(paste(path, "eqtl_analysis_misclass10p_run1.Gene", 
                                  j, ".assoc.linear.adjusted", sep = ""))
  length <- length(assoc_table[,9]) > 1
  if (length){
    for (i in 2:length(assoc_table[,9])){
      FDR <- (as.numeric(as.character(assoc_table[i,9])))
      check <- FDR < 0.05
      if (check){
        SNP_list_misclass10p[[as.character(assoc_table[i,2])]][["gene"]] <- append(SNP_list_misclass10p[[as.character(assoc_table[i,2])]][["gene"]], list(paste("Gene", j, sep = "")))
        SNP_list_misclass10p[[as.character(assoc_table[i,2])]][["pvalue"]] <- append(SNP_list_misclass10p[[as.character(assoc_table[i,2])]][["pvalue"]], list(FDR))
      }
    }
  }
}

# summary_snps_c1 <- as.data.frame(summary(SNP_list_ct1))
# length_list <- list()
# for (i in 1:length(SNP_list_ct1)){
#   length_list <- append(length_list, list(length(SNP_list_ct1[[i]][["gene"]])))
# }
# summary_snps_c1 <- summary_snps_c1[1:25575,]
# summary_snps_c1[['genes']] <- length_list
# summary_snps_c1 <- cbind(summary_snps_c1['Var1'], summary_snps_c1['genes'])
# length_snps_c1 <- as.data.frame(summary_snps_c1, colnames = c('SNP', 'siginificant genes'))
# length_snps_c1 <- as.data.frame(lapply(length_snps_c1, unlist))
# ordered_length_snps_c1 <- length_snps_c1[order(-length_snps_c1$genes),]
# 
# datatable(summary_snps_c1)
# 
# summary_snps_c2 <- as.data.frame(summary(SNP_list_ct2))
# length_list <- list()
# for (i in 1:length(SNP_list_ct2)){
#   length_list <- append(length_list, list(length(SNP_list_ct2[[i]][["gene"]])))
# }
# summary_snps_c2 <- summary_snps_c2[1:25497,]
# summary_snps_c2[['genes']] <- length_list
# summary_snps_c2 <- cbind(summary_snps_c2['Var1'], summary_snps_c2['genes'])
# length_snps_c2 <- as.data.frame(summary_snps_c2, colnames = c('SNP', 'siginificant genes'))
# length_snps_c2 <- as.data.frame(lapply(length_snps_c2, unlist))
# ordered_length_snps_c2 <- length_snps_c2[order(-length_snps_c2$genes),]

# summary_snps_c1_batches <- as.data.frame(summary(SNP_list_ct1_batches))
# length_list <- list()
# for (i in 1:length(SNP_list_ct1_batches)){
#   length_list <- append(length_list, list(length(SNP_list_ct1_batches[[i]][["gene"]])))
# }
# summary_snps_c1_batches <- summary_snps_c1_batches[1:2837,]
# summary_snps_c1_batches[['genes']] <- length_list
# summary_snps_c1_batches <- cbind(summary_snps_c1_batches['Var1'], summary_snps_c1_batches['genes'])
# length_snps_c1_batches <- as.data.frame(summary_snps_c1_batches, colnames = c('SNP', 'siginificant genes'))
# length_snps_c1_batches <- as.data.frame(lapply(length_snps_c1_batches, unlist))
# ordered_length_snps_c1_batches <- length_snps_c1_batches[order(-length_snps_c1_batches$genes),]
# 
# summary_snps_c2_batches <- as.data.frame(summary(SNP_list_ct2_batches))
# length_list <- list()
# for (i in 1:length(SNP_list_ct2_batches)){
#   length_list <- append(length_list, list(length(SNP_list_ct2_batches[[i]][["gene"]])))
# }
# summary_snps_c2_batches <- summary_snps_c2_batches[1:3015,]
# summary_snps_c2_batches[['genes']] <- length_list
# summary_snps_c2_batches <- cbind(summary_snps_c2_batches['Var1'], summary_snps_c2_batches['genes'])
# length_snps_c2_batches <- as.data.frame(summary_snps_c2_batches, colnames = c('SNP', 'siginificant genes'))
# length_snps_c2_batches <- as.data.frame(lapply(length_snps_c2_batches, unlist))
# ordered_length_snps_c2_batches <- length_snps_c2_batches[order(-length_snps_c2_batches$genes),]

summary_snps_misclass1p <- as.data.frame(summary(SNP_list_misclass1p))
length_list <- list()
for (i in 1:length(SNP_list_misclass1p)){
  length_list <- append(length_list, list(length(SNP_list_misclass1p[[i]][["gene"]])))
}
summary_snps_misclass1p <- summary_snps_misclass1p[1:2537,]
summary_snps_misclass1p[['genes']] <- length_list
summary_snps_misclass1p <- cbind(summary_snps_misclass1p['Var1'], summary_snps_misclass1p['genes'])
length_snps_misclass1p <- as.data.frame(summary_snps_misclass1p, colnames = c('SNP', 'siginificant genes'))
length_snps_misclass1p <- as.data.frame(lapply(length_snps_misclass1p, unlist))
ordered_length_snps_misclass1p <- length_snps_misclass1p[order(-length_snps_misclass1p$genes),]

summary_snps_misclass5p <- as.data.frame(summary(SNP_list_misclass5p))
length_list <- list()
for (i in 1:length(SNP_list_misclass5p)){
  length_list <- append(length_list, list(length(SNP_list_misclass5p[[i]][["gene"]])))
}
summary_snps_misclass5p <- summary_snps_misclass5p[1:2823,]
summary_snps_misclass5p[['genes']] <- length_list
summary_snps_misclass5p <- cbind(summary_snps_misclass5p['Var1'], summary_snps_misclass5p['genes'])
length_snps_misclass5p <- as.data.frame(summary_snps_misclass5p, colnames = c('SNP', 'siginificant genes'))
length_snps_misclass5p <- as.data.frame(lapply(length_snps_misclass5p, unlist))
ordered_length_snps_misclass5p <- length_snps_misclass5p[order(-length_snps_misclass5p$genes),]

summary_snps_misclass10p <- as.data.frame(summary(SNP_list_misclass10p))
length_list <- list()
for (i in 1:length(SNP_list_misclass10p)){
  length_list <- append(length_list, list(length(SNP_list_misclass10p[[i]][["gene"]])))
}
summary_snps_misclass10p <- summary_snps_misclass10p[1:2530,]
summary_snps_misclass10p[['genes']] <- length_list
summary_snps_misclass10p <- cbind(summary_snps_misclass10p['Var1'], summary_snps_misclass10p['genes'])
length_snps_misclass10p <- as.data.frame(summary_snps_misclass10p, colnames = c('SNP', 'siginificant genes'))
length_snps_misclass10p <- as.data.frame(lapply(length_snps_misclass10p, unlist))
ordered_length_snps_misclass10p <- length_snps_misclass10p[order(-length_snps_misclass10p$genes),]
