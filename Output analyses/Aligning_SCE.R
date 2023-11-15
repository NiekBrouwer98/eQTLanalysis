library(harmony)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(splatter)
library(ggplot2)
library(cowplot)

#Aligning two SingleCell Experiments
for (i in 1:length(x=sim_count_table_list)){
  A_sim_name <- paste("A_sim_", i, sep = "")
  seurat_sim <- CreateSeuratObject(sim_count_table_list[[i]])
  assign(A_sim_name, seurat_sim)
}

A_sim_list = list(A_sim_1,A_sim_2, A_sim_3, A_sim_4, A_sim_5,
                  A_sim_6, A_sim_7, A_sim_8, A_sim_9, A_sim_10)

batches = FALSE
if (batches){
  for (i in 1:5){
    A_sim_list[[i]]@meta.data[, "celltype"] <- "cell type 1"
    ID_name <- paste("ID", i, sep = "")
    A_sim_list[[i]]@meta.data[, "simulation"] <- ID_name
  }
  for (i in 6:10){
    A_sim_list[[i]]@meta.data[, "celltype"] <- "cell type 2"
    ID_name <- paste("ID", i-5, sep = "")
    A_sim_list[[i]]@meta.data[, "simulation"] <- ID_name
  }
  
}

if (batches==FALSE){
  for (i in 1:length(x=A_sim_list)){
    simulation_name <- paste("ct", i, sep = "")
    A_sim_list[[i]]@meta.data[, "simulation"] <- simulation_name
  }
  
  for (i in 1:length(x=A_sim_list)){
    A_sim_list[[i]]@meta.data[, "celltype"] <- groups_list[[i]]
  }
  
}

for (i in 1:length(x= A_sim_list)) {
  A_sim_list[[i]] <- NormalizeData(object = A_sim_list[[i]], verbose = FALSE)
  A_sim_list[[i]] <- FindVariableFeatures(object = A_sim_list[[i]], 
                                          selection.method = "vst", nfeatures = 1500, verbose = FALSE)
}
  
#Anchor Alignment
A_sim.anchors <- FindIntegrationAnchors(object.list = A_sim_list, anchor.features = 1500, dims = 1:30)
A_sim.integrated <- IntegrateData(anchorset = A_sim.anchors, dims = 1:30)


# #FastMNN alignment
# A_sim.integrated <- RunFastMNN(object.list = A_sim_list, features=700)
# A_sim.integrated <- RunUMAP(object = A_sim.integrated, reduction = "mnn", dims = 1:50)                            
# p1 <- DimPlot(object = A_sim.integrated, group.by = "simulation") + scale_color_npg()
# p2 <- DimPlot(object = A_sim.integrated, group.by = "celltype", label = FALSE, repel = TRUE)  + scale_color_npg()
# 
# plot_grid(p1, p2)



# #Harmony alignment
# sim_1and2 <- merge(A_sim_list[[1]], A_sim_list[[2]])
# sim_3and4 <- merge(A_sim_list[[3]], A_sim_list[[4]])
# sim_1234 <- merge(sim_1and2, sim_3and4)
# sim_12345 <- merge(sim_1234, A_sim_list[[5]])
# A_sim.integrated <- NormalizeData(sim_12345) %>% FindVariableFeatures(features=400) %>% ScaleData() %>% RunPCA(verbose = FALSE)
# test_p1 <- DimPlot(object = A_sim.integrated, reduction = "pca", pt.size = .1, group.by = "simulation")
# test_p2 <- VlnPlot(object = A_sim.integrated, features = "PC_1", group.by = "simulation", pt.size = .1)
# plot_grid(test_p1, test_p2)
# A_sim.integrated <- RunHarmony(A_sim.integrated, group.by.vars = "simulation", plot_convergence = TRUE)
# A_sim.integrated <- RunUMAP(object = A_sim.integrated, reduction = "harmony", dims = 1:50)
# p1 <- DimPlot(object = A_sim.integrated, group.by = "simulation") + scale_color_npg()
# p2 <- DimPlot(object = A_sim.integrated, group.by = "celltype", label = FALSE, repel = TRUE)  + scale_color_npg()
# 
# plot_grid(p1, p2)

  
# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
DefaultAssay(object = A_sim.integrated) <- "integrated"


# Run the standard workflow for visualization and clustering
A_sim.integrated <- ScaleData(object = A_sim.integrated, verbose = FALSE)
A_sim.integrated <- RunPCA(object = A_sim.integrated, npcs = 50, verbose = FALSE)
A_sim.integrated <- RunUMAP(object = A_sim.integrated, reduction = "pca", 
                               dims = 1:30)
A_sim.integrated <- FindNeighbors(A_sim.integrated, reduction = "pca", dims = 1:30)
A_sim.integrated <- FindClusters(A_sim.integrated, resolution = 0.1)

# Visualisation
# library(ggsci)
# 
# p1 <- DimPlot(object = A_sim.integrated, reduction = "umap", group.by = "simulation") + scale_color_npg()
# p2 <- DimPlot(object = A_sim.integrated, reduction = "umap", group.by = "celltype",
#               label = FALSE, repel = TRUE)  + scale_color_npg()
# 
# plot_grid(p1, p2)

# Classification
# A_sim.integrated@meta.data[["integrated_snn_res.0.3"]]


A_sim.integrated@meta.data[,"simulation_number"] <- "initialize"
for (i in 1:length(x=A_sim.integrated@meta.data[,"celltype"])){
  if (i < 2001){
    if (A_sim.integrated@meta.data[,"celltype"][i] == "Group1"){
      A_sim.integrated@meta.data[,"simulation_number"][i] <- "ID1.t1"
    }
    if (A_sim.integrated@meta.data[,"celltype"][i] == "Group2"){
      A_sim.integrated@meta.data[,"simulation_number"][i] <- "ID1.t2"
    }
  }
  if (i > 2000){
    if (A_sim.integrated@meta.data[,"celltype"][i] == "Group1"){
      A_sim.integrated@meta.data[,"simulation_number"][i] <- "ID2.t1"
    }
    if (A_sim.integrated@meta.data[,"celltype"][i] == "Group2"){
      A_sim.integrated@meta.data[,"simulation_number"][i] <- "ID2.t2"
    }
  }    
}

A_sim.integrated@meta.data[,"prediction_number"] <- "initialize"
for (i in 1:length(x=A_sim.integrated@meta.data[,"integrated_snn_res.0.1"])){
  if (i < 2001){
    if (A_sim.integrated@meta.data[,"integrated_snn_res.0.1"][i] == "0"){
      A_sim.integrated@meta.data[,"prediction_number"][i] <- "ID1.t1"
    }
    if (A_sim.integrated@meta.data[,"integrated_snn_res.0.1"][i] == "1"){
      A_sim.integrated@meta.data[,"prediction_number"][i] <- "ID1.t2"
    }
  }
  if (i > 2000){
    if (A_sim.integrated@meta.data[,"integrated_snn_res.0.1"][i] == "0"){
      A_sim.integrated@meta.data[,"prediction_number"][i] <- "ID2.t1"
    }
    if (A_sim.integrated@meta.data[,"integrated_snn_res.0.1"][i] == "1"){
      A_sim.integrated@meta.data[,"prediction_number"][i] <- "ID2.t2"
    }
  }    
}


A_sim.integrated$prediction.match <- A_sim.integrated@meta.data[,"simulation_number"] == A_sim.integrated@meta.data[, "prediction_number"]
table(A_sim.integrated$prediction.match)
table(real = A_sim.integrated@meta.data[,"simulation_number"], pred = A_sim.integrated@meta.data[,"prediction_number"])
