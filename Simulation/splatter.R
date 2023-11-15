library(splatter)
library(scater)
library(fitdistrplus)
library(ggplot2)
library(cowplot)
library(ggsci)
library(ggpubr)

celltype1 <- rowMeans(Y_nl)
celltype2 <- rowMeans(Y_celltype2_nl)

fit.gamma_ct1 <- fitdist(celltype1, distr = "gamma", method = 'mle')
fit.gamma_ct2 <- fitdist(celltype2, distr = "gamma", method = 'mle')
shape_ct1 <- fit.gamma_ct1$estimate['shape']
rate_ct1 <- fit.gamma_ct1$estimate['rate']
shape_ct2 <- fit.gamma_ct2$estimate['shape']
rate_ct2 <- fit.gamma_ct2$estimate['rate']

#Determine gamma shape and rate per ID
for (i in 1:10){
  shape <- paste("shape", i, sep = "")
  rate <- paste("rate", i, sep = "")
  
  fit.gamma <- fitdist(Y_nl[i,], distr = "gamma", method = 'mle')
  summary(fit.gamma)
  assign(shape, fit.gamma$estimate['shape'])
  assign(rate, fit.gamma$estimate['rate'])
  # plot(fit.gamma, demp=TRUE)
  
}

#Setting the parameters
params <- newSplatParams()
params
params_easy <- setParams(params, nGenes = 1500, batchCells = c(2000),
                        group.prob = c(0.75, 0.25),
                          de.prob = 0.3, de.downProb = c(0.5,0.1),
                          de.facLoc = 0.3, de.facScale = 0.01)

params_int <- setParams(params, nGenes = 1500, batchCells = c(2000),
                         group.prob = c(0.75, 0.25),
                         de.prob = 0.1, de.downProb = c(0.5, 0.1),
                        de.facLoc = 0.3, de.facScale = 0.01)

params_hard <- setParams(params, nGenes = 1500, batchCells = c(2000),
                         de.prob = 0.01,
                         group.prob = c(0.75, 0.25),
                         de.prob = 0.01, de.downProb = c(0.5, 0.1),
                         de.facLoc = 0.3, de.facScale = 0.01)

params_batches <- setParams(params, nGenes = 1500,
                            batch.facLoc = 0.01, batch.facScale = 0.01)
                            # group.prob = c(0.75, 0.25),
                            # de.prob = 0.3, de.downProb = c(0.5, 0.1),
                            # de.facLoc = 0.3, de.facScale = 0.01)


#Explore how this looks
visualise_scenarios = FALSE
if (visualise_scenarios) {
  sim_easy = splatSimulate(params_easy, method = "groups")
  sim_easy_norm = normalizeSCE(sim_easy)
  sim_int = splatSimulate(params_int, method = "groups")
  sim_int_norm = normalizeSCE(sim_int)
  sim_hard = splatSimulate(params_hard, method = "groups")
  sim_hard_norm = normalizeSCE(sim_hard)
  
  p1 <- plotTSNE(sim_easy_norm, colour_by ="Group")
  p2 <- plotTSNE(sim_int_norm, colour_by ="Group")
  p3 <- plotTSNE(sim_hard_norm, colour_by ="Group")
  tsne1 <- ggplot(p1[["data"]], aes(X, Y, color= p1[["data"]][["colour_by"]])) + geom_point() + scale_color_npg() + 
    labs(color = "Groups") + xlab("Dimension 1") + ylab("Dimension 2")
  tsne2 <- ggplot(p2[["data"]], aes(X, Y, color= p2[["data"]][["colour_by"]])) + geom_point() + scale_color_npg() + 
    labs(color = "Groups") + xlab("Dimension 1") + ylab("Dimension 2")
  tsne3 <- ggplot(p3[["data"]], aes(X, Y, color= p3[["data"]][["colour_by"]])) + geom_point() + scale_color_npg() + 
    labs(color = "Groups") + xlab("Dimension 1") + ylab("Dimension 2")
  
  ggarrange(tsne1, tsne2, tsne3, ncol = 3, nrow = 1, common.legend = T)
}

#generate for each individual
params1 <- setParams(params_easy, update = list(mean.shape = shape1, mean.rate= rate1))
params2 <- setParams(params_easy, update = list(mean.shape = shape2, mean.rate= rate2))
params3 <- setParams(params_easy, update = list(mean.shape = shape3, mean.rate= rate3))
params4 <- setParams(params_easy, update = list(mean.shape = shape4, mean.rate= rate4))
params5 <- setParams(params_easy, update = list(mean.shape = shape5, mean.rate= rate5))
params6 <- setParams(params_easy, update = list(mean.shape = shape6, mean.rate= rate6))
params7 <- setParams(params_easy, update = list(mean.shape = shape7, mean.rate= rate7))
params8 <- setParams(params_easy, update = list(mean.shape = shape8, mean.rate= rate8))
params9 <- setParams(params_easy, update = list(mean.shape = shape9, mean.rate= rate9))
params10 <- setParams(params_easy, update = list(mean.shape = shape10, mean.rate= rate10))
params_shape <- setParams(params_easy, update = list(mean.shape = shape1, mean.rate= rate2))
params_rate <- setParams(params_easy, update = list(mean.shape = shape2, mean.rate= rate1))
params_0.01 <- setParams(params_easy, update = list(mean.shape = 0.01+shape1, mean.rate= rate1))
params_0.1 <- setParams(params_easy, update = list(mean.shape = 0.1+shape1, mean.rate= rate1))
params_0.5 <- setParams(params_easy, update = list(mean.shape = 0.5+shape1, mean.rate= rate1))
params_1 <- setParams(params_easy, update = list(mean.shape = 1+shape1, mean.rate= rate1))

fifty_list_c1 <- c()
fifty_list_c2 <- c()
for (i in 1:50){
  fifty_list_c1 <- append(fifty_list_c1, 300)
  fifty_list_c2 <- append(fifty_list_c2,100)
}


params_ct1 <- setParams(params_batches, update = list(mean.shape = shape_ct1, mean.rate= rate_ct1, batchCells = fifty_list_c1))
params_ct2 <- setParams(params_batches, update = list(mean.shape = shape_ct2, mean.rate= rate_ct2, batchCells = fifty_list_c2))

# params_list <- list(params1, params_0.01, params_0.1, params_0.5, params_1)
# params_list <- list(params1, params2, params3, params4, params5, params6, params7, params8, params9, params10)
params_list <- list(params_ct1, params_ct2)

#Simulating the experiments
sim_list <- list()
for (i in 1:length(x=params_list)){
  sim_list[[i]] <- splatSimulate(params_list[[i]], method = "group")
  }

#Generate count table list
sim_count_table_list <- list()
for (i in 1:length(x=sim_list)){
  sim_count_table_list[[i]] <- counts(sim_list[[i]])
}

#Visulaise the first two sims seperately
sim1_norm <- normalizeSCE(sim_list[[1]])
sim2_norm <- normalizeSCE(sim_list[[2]])
color_scale <- scale_colour_manual(name = "Batch", values = scale_color_npg())
tsne_sim1 <- plotTSNE(sim1_norm, colour_by= "Batch")
ggplot(tsne_sim1[["data"]], aes(X, Y, color= tsne_sim1[["data"]][["colour_by"]])) + geom_point(alpha=0.5) + scale_color_npg() +
  labs(title = 't-sne of cell type 1', color = "Individuals") + xlab("Dimension 1") + ylab("Dimension 2")
tsne_sim2 <- plotTSNE(sim2_norm, colour_by= "Batch")
ggplot(tsne_sim2[["data"]], aes(X, Y, color= tsne_sim2[["data"]][["colour_by"]])) + geom_point(alpha=0.5) + scale_color_npg() +
  labs(title = 't-sne of cell type 2', color = "Individuals") + xlab("Dimension 1") + ylab("Dimension 2")

library(scMerge)

rm(A_sim.anchors)

SC_combine <- sce_cbind(list(sim1_norm, sim2_norm), batch_names = c("cell type 1","cell type 2"))
tsne_sc_combine <- plotTSNE(SC_combine, colour_by= "batch")
ggplot(tsne_sc_combine[["data"]], aes(X, Y, color= tsne_sc_combine[["data"]][["colour_by"]])) + geom_point(alpha=0.03) + scale_color_npg() +
  labs(title = 't-sne of merged cell types', color = "") + xlab("Dimension 1") + ylab("Dimension 2")

batches = TRUE
if (batches){
  start_c1 <- seq(1, by = 300, length = ncol(sim_count_table_list[[1]]) / 300)
  start_c2 <- seq(1, by = 100, length = ncol(sim_count_table_list[[2]]) / 100)
  IDs_list_c1 <- lapply(start_c1, function(i, df) df[,i:(i+299)], df = sim_count_table_list[[1]])
  IDs_list_c2 <- lapply(start_c2, function(i, df) df[,i:(i+99)], df = sim_count_table_list[[2]])
  #Set column names
  for (i in 1:length(x=IDs_list_c1)){
    for (j in 1:300){
      name <- paste("ID", i, "_cell_type1_", j, sep = "")
      colnames(IDs_list_c1[[i]])[j] <- name
    }
  }
  for (i in 1:length(x=IDs_list_c2)){
    for (j in 1:100){
      name <- paste("ID", i, "_cell_type2_", j, sep = "")
      colnames(IDs_list_c2[[i]])[j] <- name
    }
  }
}

if (batches==FALSE) {
  #Generate group label list
  groups_list <- list()
  for (i in 1:length(x=sim_list)){
    groups_list[[i]] <- colData(sim_list[[i]])$Group
  }
  
  #Set column names
  for (i in 1:length(x=sim_list)){
    for (j in 1:2000){
      name <- paste("ID", i, "_cell", j, sep = "")
      colnames(sim_count_table_list[[i]])[j] <- name
    }
  }
}

if (exists('sim_list')& exists('groups_list')& exists('sim_count_table_list')){
  print("The result are 3 lists: sim_list, groups_list (with all the celltypes),
        sim_count_table_list (all count tables with named columns)")
}
