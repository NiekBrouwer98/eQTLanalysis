is.contained=function(vec1,vec2){
  x=vector(length = length(vec1))
  for (i in 1:length(vec1)) {
    x[i] = vec1[i] %in% vec2
    if(length(which(vec1[i] %in% vec2)) == 0) vec2 else 
      vec2=vec2[-match(vec1[i], vec2)]
  }
  y=all(x==T)
  return(y)
}

snps_list <- list()
for (i in 1:nrow(ordered_length_snps_misclass1p['Var1'])){
  name <- ordered_length_snps_misclass1p[i,'Var1']
  snps <- genotypes50[[as.character(name)]]
  snps_list <- append(snps_list, list(snps))
}
ordered_length_snps_misclass1p$snps <- snps_list

filtered_ordered_length_snps_misclass1p <- data.frame()
for (i in 1:nrow(ordered_length_snps_misclass1p)){
  check <- is.contained(c(0,1,2), ordered_length_snps_misclass1p[i,'snps'][[1]])
  if (check==TRUE){
    filtered_ordered_length_snps_misclass1p <- rbind(filtered_ordered_length_snps_misclass1p, ordered_length_snps_misclass1p[i,])
  }
}

snps_list <- list()
for (i in 1:nrow(ordered_length_snps_misclass5p['Var1'])){
  name <- ordered_length_snps_misclass5p[i,'Var1']
  snps <- genotypes50[[as.character(name)]]
  snps_list <- append(snps_list, list(snps))
}
ordered_length_snps_misclass5p$snps <- snps_list

filtered_ordered_length_snps_misclass5p <- data.frame()
for (i in 1:nrow(ordered_length_snps_misclass5p)){
  check <- is.contained(c(0,1,2), ordered_length_snps_misclass5p[i,'snps'][[1]])
  if (check==TRUE){
    filtered_ordered_length_snps_misclass5p <- rbind(filtered_ordered_length_snps_misclass5p, ordered_length_snps_misclass5p[i,])
  }
}

snps_list <- list()
for (i in 1:nrow(ordered_length_snps_misclass10p['Var1'])){
  name <- ordered_length_snps_misclass10p[i,'Var1']
  snps <- genotypes50[[as.character(name)]]
  snps_list <- append(snps_list, list(snps))
}
ordered_length_snps_misclass10p$snps <- snps_list

filtered_ordered_length_snps_misclass10p <- data.frame()
for (i in 1:nrow(ordered_length_snps_misclass10p)){
  check <- is.contained(c(0,1,2), ordered_length_snps_misclass10p[i,'snps'][[1]])
  if (check==TRUE){
    filtered_ordered_length_snps_misclass10p <- rbind(filtered_ordered_length_snps_misclass10p, ordered_length_snps_misclass10p[i,])
  }
}






snps_list <- list()
for (i in 1:nrow(ordered_length_snps_c1['Var1'])){
  name <- ordered_length_snps_c1[i,'Var1']
  snps <- genotypes[[as.character(name)]]
  snps_list <- append(snps_list, list(snps))
}
ordered_length_snps_c1$snps <- snps_list

filtered_ordered_length_snps_c1 <- data.frame()
for (i in 1:nrow(ordered_length_snps_c1)){
  check <- is.contained(c(0,1,2), ordered_length_snps_c1[i,'snps'][[1]])
  if (check==TRUE){
    filtered_ordered_length_snps_c1 <- rbind(filtered_ordered_length_snps_c1, ordered_length_snps_c1[i,])
  }
}

snps_list <- list()
for (i in 1:nrow(ordered_length_snps_c2['Var1'])){
  name <- ordered_length_snps_c2[i,'Var1']
  snps <- genotypes[[as.character(name)]]
  snps_list <- append(snps_list, list(snps))
}

ordered_length_snps_c2$snps <- snps_list

filtered_ordered_length_snps_c2 <- data.frame()
for (i in 1:nrow(ordered_length_snps_c2)){
  check <- is.contained(c(0,1,2), ordered_length_snps_c2[i,'snps'][[1]])
  if (check==TRUE){
    filtered_ordered_length_snps_c2 <- rbind(filtered_ordered_length_snps_c2, ordered_length_snps_c2[i,])
  }
}

snps_list <- list()
for (i in 1:nrow(ordered_length_snps_c1_batches['Var1'])){
  name <- ordered_length_snps_c1_batches[i,'Var1']
  snps <- genotypes50[[as.character(name)]]
  snps_list <- append(snps_list, list(snps))
}

ordered_length_snps_c1_batches$snps <- snps_list

filtered_ordered_length_snps_c1_batches <- data.frame()
for (i in 1:nrow(ordered_length_snps_c1_batches)){
  check <- is.contained(c(0,1,2), ordered_length_snps_c1_batches[i,'snps'][[1]])
  if (check==TRUE){
    filtered_ordered_length_snps_c1_batches <- rbind(filtered_ordered_length_snps_c1_batches, ordered_length_snps_c1_batches[i,])
  }
}

snps_list <- list()
for (i in 1:nrow(ordered_length_snps_c2_batches['Var1'])){
  name <- ordered_length_snps_c2_batches[i,'Var1']
  snps <- genotypes50[[as.character(name)]]
  snps_list <- append(snps_list, list(snps))
}

ordered_length_snps_c2_batches$snps <- snps_list

filtered_ordered_length_snps_c2_batches <- data.frame()
for (i in 1:nrow(ordered_length_snps_c2_batches)){
  check <- is.contained(c(0,1,2), ordered_length_snps_c2_batches[i,'snps'][[1]])
  if (check==TRUE){
    filtered_ordered_length_snps_c2_batches <- rbind(filtered_ordered_length_snps_c2_batches, ordered_length_snps_c2_batches[i,])
  }
}
