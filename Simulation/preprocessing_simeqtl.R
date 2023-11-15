library(Seurat)
library(tibble)

output_path <- 'C:/Users/niekb/Documents/BEP datasets/simulated_data'

#genotypes
path <- 'C:/Users/niekb/plink_workspace/genotype_sim'
plink_genotypes <- readStandardGenotypes(50, path, format = "plink")

#count data
#Split celltypes
# celltype1 <- subset(A_sim.integrated, subset = celltype < "Group2")
# celltype2 <- subset(A_sim.integrated, subset = celltype > "Group1")
celltype1 <- sim_count_table_list[[1]]
celltype2 <- sim_count_table_list[[2]]

#Find IDs as named in vcf file
# names_celltype1 <- list(celltype1@assays[["RNA"]]@data@Dimnames[[2]])
# names_celltype2 <- list(celltype2@assays[["RNA"]]@data@Dimnames[[2]])
# vcfnames <- genotypes[["id_samples"]][1:10]
# vcfnames_celltype1 <- list()
# for (i in 1:length(names_celltype1[[1]])){
#   pers <- strtoi(substr(names_celltype1[[1]][i], 3,3))
#   if (substr(names_celltype1[[1]][i], 4,4) == '0'){
#    vcfnames_celltype1[[i]] <- vcfnames[[10]] 
#   }
#   else{
#     vcfnames_celltype1[[i]] <- vcfnames[[pers]]  
#   }
# }
# 
# vcfnames_celltype2 <- list()
# for (i in 1:length(names_celltype2[[1]])){
#   pers <- strtoi(substr(names_celltype2[[1]][i], 3,3))
#   if (substr(names_celltype2[[1]][i], 4,4) == '0'){
#     vcfnames_celltype2[[i]] <- vcfnames[[10]] 
#   }
#   else{
#     vcfnames_celltype2[[i]] <- vcfnames[[pers]]  
#   }
# }

strtoi(substr(sim_list[[1]]@colData@listData[["Batch"]][1],6,7))

#for batches:
vcfnames <- genotypes[["id_samples"]][1:50]
vcfnames_celltype1 <- list()
for (i in 1:15000){
  pers <- strtoi(substr(sim_list[[1]]@colData@listData[["Batch"]][i],6,7))
  vcfnames_celltype1[[i]] <- vcfnames[[pers]]
}

vcfnames_celltype2 <- list()
for (i in 1:5000){
  pers <- strtoi(substr(sim_list[[2]]@colData@listData[["Batch"]][i],6,7))
  vcfnames_celltype2[[i]] <- vcfnames[[pers]]
}

#Create matrices
count_matrix_celltype1 <- as.data.frame(celltype1@assays[["RNA"]]@counts)
count_matrix_celltype2 <- as.data.frame(celltype2@assays[["RNA"]]@counts)

# colnames(count_matrix_celltype1) <- c(vcfnames_celltype1)
# colnames(count_matrix_celltype2) <- c(vcfnames_celltype2)
# rownames(count_matrix_celltype1) <- gene_to_snp[,2]
# rownames(count_matrix_celltype2) <- gene_to_snp[,2]

gene_to_snp <- data.frame(read_name= celltype1@assays[["RNA"]]@data@Dimnames[[1]],
                          SNP_name= colnames(causalSNPs))

output_data_celltype1 <- t(count_matrix_celltype1)
output_data_celltype2 <- t(count_matrix_celltype2)

#Add FID and IID columns
# dim(output_data_celltype1)
# output_data_celltype1 <- cbind(vcfnames_celltype1, vcfnames_celltype1, output_data_celltype1)
# dim(output_data_celltype1)
# # ouptut_data_celltype1 <-cbind(chr_list1, output_data_celltype1)
# output_data_celltype2 <- cbind(vcfnames_celltype2, vcfnames_celltype2, output_data_celltype2)
# dim(output_data_celltype2)

#for batches:
data_celltype1_batches <- t(celltype1)
length(vcfnames_celltype1)
dim(data_celltype1_batches)
output_data_celltype1_batches <- cbind(vcfnames_celltype1, vcfnames_celltype1, data_celltype1_batches)
dim(output_data_celltype1_batches)

data_celltype2_batches <- t(celltype2)
length(vcfnames_celltype2)
dim(data_celltype2_batches)
output_data_celltype2_batches <- cbind(vcfnames_celltype2, vcfnames_celltype2, data_celltype2_batches)
dim(output_data_celltype2_batches)

snp_name <- as.list(gene_to_snp[,1])
gene_name <- colnames(data_celltype1_batches)
length(gene_name)

column_names <- list(as.name("FID"), as.name("IID"))
for (i in 1:length(x=gene_name)){
  column_names[[2+i]] <- as.name(gene_name[i])
}

length(column_names)
dim(output_data_celltype1_batches)
colnames(output_data_celltype1_batches) <- column_names
colnames(output_data_celltype2_batches) <- column_names
rownames(output_data_celltype1_batches) <- c()
rownames(output_data_celltype2_batches) <- c()

output_data_ct1_batches_misclass_1p <- output_data_celltype1_batches
output_data_ct1_batches_misclass_5p <- output_data_celltype1_batches
output_data_ct1_batches_misclass_10p <- output_data_celltype1_batches

random_samples_1p <- sample(1:5000, 150, replace=F)
random_samples_5p <- sample(1:5000, 750, replace=F)
random_samples_10p <- sample(1:5000, 1500, replace=F)
for (i in random_samples_1p){
  output_data_ct1_batches_misclass_1p <- rbind(output_data_ct1_batches_misclass_1p,
                                            output_data_celltype2_batches[i,])
}
for (i in random_samples_5p){
  output_data_ct1_batches_misclass_5p <- rbind(output_data_ct1_batches_misclass_5p,
                                               output_data_celltype2_batches[i,])
}
for (i in random_samples_10p){
  output_data_ct1_batches_misclass_10p <- rbind(output_data_ct1_batches_misclass_10p,
                                               output_data_celltype2_batches[i,])
}

dim(output_data_ct1_batches_misclass_1p)
dim(output_data_ct1_batches_misclass_5p)
dim(output_data_ct1_batches_misclass_10p)

#Write to txt file
write.table(output_data_celltype1_batches, file=paste(output_path, '/matrix_celltype1_batches.txt', sep = ""),sep = "\t", row.names = F, quote = F)
write.table(output_data_celltype2_batches, file=paste(output_path, '/matrix_celltype2_batches.txt', sep = ""),sep = "\t", row.names = F, quote = F)

write.table(output_data_ct1_batches_misclass_1p, file=paste(output_path, '/matrix_celltype1_batches_misclass1p.txt', sep = ""),sep = "\t", row.names = F, quote = F)
write.table(output_data_ct1_batches_misclass_5p, file=paste(output_path, '/matrix_celltype1_batches_misclass5p.txt', sep = ""),sep = "\t", row.names = F, quote = F)
write.table(output_data_ct1_batches_misclass_10p, file=paste(output_path, '/matrix_celltype1_batches_misclass10p.txt', sep = ""),sep = "\t", row.names = F, quote = F)


#Write extra name files
lapply(names_celltype2, function(x) write.table( data.frame(x),paste(output_path, '/cellnames_celltype2.csv', sep = ""), append= T, sep=',' ))
write.table(gene_to_snp, paste(output_path, '/read_to_location.csv', sep = ""), append= T, sep=',' )

first_50_genotypes <- genotypes[["genotypes"]][1:50,]

write.table(first_50_genotypes, paste(output_path,'/genotypes50.txt', sep=""), sep= '\t')
