gene_count_c1 <- 0
gene_list_c1 <- list()
for (i in 1:nrow(filtered_ordered_length_snps_c1)){
  SNP <- filtered_ordered_length_snps_c1[i,'Var1']
  genes <- SNP_list_ct1[[as.character(SNP)]][["gene"]]
  for (j in 1:length(genes)){
    gene <- genes[[j]]
    check <- (as.character(gene) %in% gene_list_c1)
    if (check==FALSE){
      gene_list_c1 <- append(gene_list_c1, as.character(gene))
      gene_count_c1 <- gene_count_c1 +1
    }
  }
}

gene_count_c2 <- 0
gene_list_c2 <- list()
for (i in 1:nrow(filtered_ordered_length_snps_c2)){
  SNP <- filtered_ordered_length_snps_c2[i,'Var1']
  genes <- SNP_list_ct2[[as.character(SNP)]][["gene"]]
  for (j in 1:length(genes)){
    gene <- genes[[j]]
    check <- (as.character(gene) %in% gene_list_c2)
    if (check==FALSE){
      gene_list_c2 <- append(gene_list_c2, as.character(gene))
      gene_count_c2 <- gene_count_c2 +1
    }
  }
}     

gene_count_c1_batches <- 0
gene_list_c1_batches <- list()
for (i in 1:nrow(filtered_ordered_length_snps_c1_batches)){
  SNP <- filtered_ordered_length_snps_c1_batches[i,'Var1']
  genes <- SNP_list_ct1_batches[[as.character(SNP)]][["gene"]]
  for (j in 1:length(genes)){
    gene <- genes[[j]]
    check <- (as.character(gene) %in% gene_list_c1_batches)
    if (check==FALSE){
      gene_list_c1_batches <- append(gene_list_c1_batches, as.character(gene))
      gene_count_c1_batches <- gene_count_c1_batches +1
    }
  }
}

gene_count_c2_batches <- 0
gene_list_c2_batches <- list()
for (i in 1:nrow(filtered_ordered_length_snps_c2_batches)){
  SNP <- filtered_ordered_length_snps_c2_batches[i,'Var1']
  genes <- SNP_list_ct2_batches[[as.character(SNP)]][["gene"]]
  for (j in 1:length(genes)){
    gene <- genes[[j]]
    check <- (as.character(gene) %in% gene_list_c2_batches)
    if (check==FALSE){
      gene_list_c2_batches <- append(gene_list_c2_batches, as.character(gene))
      gene_count_c2_batches <- gene_count_c2_batches +1
    }
  }
}   

length(intersect(filtered_ordered_length_snps_c1['Var1'][[1]], filtered_ordered_length_snps_c2['Var1'][[1]]))
(intersect(filtered_ordered_length_snps_c1_batches['Var1'][[1]], filtered_ordered_length_snps_c2_batches['Var1'][[1]]))
length(intersect(gene_list_c1, gene_list_c2))
(intersect(gene_list_c1_batches, gene_list_c2_batches))[[1]]

gene_count_misclass1p <- 0
gene_list_misclass1p <- list()
for (i in 1:nrow(filtered_ordered_length_snps_misclass1p)){
  SNP <- filtered_ordered_length_snps_misclass1p[i,'Var1']
  genes <- SNP_list_misclass1p[[as.character(SNP)]][["gene"]]
  for (j in 1:length(genes)){
    gene <- genes[[j]]
    check <- (as.character(gene) %in% gene_list_misclass1p)
    if (check==FALSE){
      gene_list_misclass1p <- append(gene_list_misclass1p, as.character(gene))
      gene_count_misclass1p <- gene_count_misclass1p +1
    }
  }
} 

gene_count_misclass5p <- 0
gene_list_misclass5p <- list()
for (i in 1:nrow(filtered_ordered_length_snps_misclass5p)){
  SNP <- filtered_ordered_length_snps_misclass5p[i,'Var1']
  genes <- SNP_list_misclass5p[[as.character(SNP)]][["gene"]]
  for (j in 1:length(genes)){
    gene <- genes[[j]]
    check <- (as.character(gene) %in% gene_list_misclass5p)
    if (check==FALSE){
      gene_list_misclass5p <- append(gene_list_misclass5p, as.character(gene))
      gene_count_misclass5p <- gene_count_misclass5p +1
    }
  }
} 

gene_count_misclass10p <- 0
gene_list_misclass10p <- list()
for (i in 1:nrow(filtered_ordered_length_snps_misclass10p)){
  SNP <- filtered_ordered_length_snps_misclass10p[i,'Var1']
  genes <- SNP_list_misclass10p[[as.character(SNP)]][["gene"]]
  for (j in 1:length(genes)){
    gene <- genes[[j]]
    check <- (as.character(gene) %in% gene_list_misclass10p)
    if (check==FALSE){
      gene_list_misclass10p <- append(gene_list_misclass10p, as.character(gene))
      gene_count_misclass10p <- gene_count_misclass10p +1
    }
  }
} 

length(intersect(filtered_ordered_length_snps_c1_batches['Var1'][[1]], filtered_ordered_length_snps_misclass1p['Var1'][[1]]))
length(intersect(gene_list_c1_batches, gene_list_misclass1p))

length(intersect(filtered_ordered_length_snps_c1_batches['Var1'][[1]], filtered_ordered_length_snps_misclass5p['Var1'][[1]]))
length(intersect(gene_list_c1_batches, gene_list_misclass5p))

length(intersect(filtered_ordered_length_snps_c1_batches['Var1'][[1]], filtered_ordered_length_snps_misclass10p['Var1'][[1]]))
length(intersect(gene_list_c1_batches, gene_list_misclass10p))

