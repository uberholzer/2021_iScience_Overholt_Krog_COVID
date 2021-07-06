library(Seurat)
library(ggplot2)
library(future)
library(dplyr)
library(harmony)
set.seed(10097)

print(packageVersion("Seurat"))
update.packages("Seurat")
print(packageVersion("Seurat"))


# # Read in coarse annotations
# immune.combined <- readRDS("main_coarse_labeled.rds")
# 
# # Set clustering paramters
# PCTthresh <- 0.1 # Set threshold here for determining the number of PCs to use
# n_beyond <- 5
# 
# 
# ########################## PCB CLUSTERING ###########################
# # Include ''Plasma cell', 'B cell', and 'Proliferating' annotations
# 
# Idents(object=immune.combined) <- "celltype_coarse"
# PCB = subset(immune.combined,idents = c("Plasma cell","B cell","Proliferating")) # Get plasma cell, B cell, and proliferating clusters
# DefaultAssay(PCB) <- "RNA"
# 
# # Conduct test to choose the number of PCs to include
# pca.test<- PCB %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))
# 
# var.genes <- pca.test@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', pca.test@assays$RNA@var.features)]
# pca.test <- RunPCA(pca.test, features = var.genes, npcs = 100, verbose = FALSE) # Start by calculating 100 PCs
# 
# # Choose number of PCs
# pca.test <- ProjectDim(object = pca.test)
# pdf(file="PCB1+5_elbowplot.pdf")
# ElbowPlot(object = pca.test,ndims = 100)
# dev.off()
# 
# pct <- pca.test[["pca"]]@stdev / sum(pca.test[["pca"]]@stdev) * 100
# nPCs <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > PCTthresh), decreasing = T)[1] + 1 + n_beyond
# 
# cat("/nNumber of PCs used for PCB: ")
# print(nPCs)
# 
# # Run the actual PCA using the number of PCs calculated
# PCB <- PCB %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))
# 
# removed.genes <- PCB@assays$RNA@var.features[grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', PCB@assays$RNA@var.features)]
# print(paste("Removed genes: ",removed.genes))
# 
# var.genes <- PCB@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', PCB@assays$RNA@var.features)]
# PCB <- RunPCA(PCB, features = var.genes, npcs = nPCs, verbose = FALSE)
# 
# print("Scale and PCA successful!")
# 
# # Run Harmony
# PCB <- PCB %>% RunHarmony("orig.ident", plot_convergence = FALSE) # make sure to specify what variable(s) to regress out here
# 
# # Clustering
# PCB <- PCB %>%
#   RunUMAP(reduction = "harmony", dims = 1:nPCs) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:nPCs) %>%
#   FindClusters(resolution = 0.5) %>%
#   identity()
# 
# print("\nTotal number of cells: ")
# print(ncol(PCB@assays[["RNA"]]))
# 
# # Visualization
# pdf(file="PCBPro1+5_clusters.pdf")
# DimPlot(PCB, reduction = "umap", label = TRUE, pt.size = .1)
# dev.off()
# 
# pdf(file="PCBPro1+5_donors.pdf",width=18,height=6)
# DimPlot(object = PCB, reduction = "umap", pt.size = .1, group.by = "orig.ident",label=FALSE)
# dev.off()
# 
# pdf(file="PCBPro1+5_study.pdf",width=6,height=6)
# DimPlot(object = PCB, reduction = "umap", pt.size = .1, group.by = "study",label=FALSE)
# dev.off()
# 
# saveRDS(PCB,"PCBPro1+5.rds")


########################## MARKER DOTPLOTS AND HEATMAPS ###########################


PCB <- readRDS("PCBPro1+5.rds")

Idents(object=PCB) <- "seurat_clusters"
DefaultAssay(PCB) <- "RNA"

# Dotplots for cell type annotation
features=rev(c("TPPP3","SCGB1A1", # Epithelial
               "TYMS","MKI67", # Proliferating
               "IGHA1","IGHG4", # Plasma cells
               "IGHM","IGHD","MS4A1", # B cells
               "NKG7","KLRD1", # NK/effector
               "CD8A","CD4","CD3D", # T cells
               "LILRA4","IL3RA","CD1C", # Dendritic cells
               "FCGR3A","CD14","CD68","CD16","VCAN", # Monocyte/macrophages
               "PPBP", # Platelets
               "HBB")) # RBCs

p4 <- DotPlot(PCB, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "PCBPro_dotplot_simple.pdf",height=10,width=15)
plot(p4)
p4
dev.off()

features=rev(c("EMP2","TPPP3","KRT18","KRT15","SCGB1A1","TPSB2","TPSAB1","TYMS","MKI67","IGHA1","IGHG4","IGHM",
               "IGHD","CD22","CD19","MS4A1","CD79A","CXCR6","GZMK","KLRB1","EOMES","KLRF1","GZMB","GNLY","KLRD1",
               "CD69","FAS","FOXP3","CTLA4","IL2RA","CD44","IL7R","IL2RG","TCF7","CCR7","SELL","CD8A","CD4","CD3D","LY6GD",
               "CXCR2","FCGR3B","HCAR3","PI3","ELANE", "LILRA4","IL3RA","THBD","CD1C","CLEC9A","XCR1","FABP4","SPP1",
               "FCGR3A","CD14","S100A8","FCN1","CD68","PPBP","HBB"))

p4 <- DotPlot(PCB, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "PCBPro1+5_dotplot.pdf",height=10,width=36)
plot(p4)
p4
dev.off()

# Make heatmap
# First scale all genes to make a complete heatmap
# PCB <- ScaleData(PCB, features = rownames(PCB))

# # Parallelize FindMarkers step
# plan("multiprocess") # works with 4 workers
# plan()
# options(future.globals.maxSize = 10000 * 1024^9,seed=TRUE)
# 
# # find markers for every cluster compared to all remaining cells, report only the positive ones
# PCB@misc$markers <- FindAllMarkers(object = PCB, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
# write.table(PCB@misc$markers,file='PCB1+5_markers.txt',row.names = FALSE,quote = FALSE,sep = '\t')

##add average expression information
# PCB@misc$averageExpression = AverageExpression(object = PCB)
# write.table(PCB@misc$averageExpression$RNA,file='PCB-average_MAST.txt',row.names = TRUE,quote = FALSE,sep = '\t')

# hc.markers = read.delim2("PCB1+5wPro_markers.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
# hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
# 
# tt = DoHeatmap(object = PCB, features = top10$gene) + NoLegend()  +
#   scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdPu"))
# ggplot2::ggsave(file="PCB1+5wPro_heatmap.pdf",plot = tt,device = 'pdf',width = 16, height = 12, units = "in",dpi = 300)



########################## CLUSTER ANNOTATION ###########################

Idents(object=PCB) <- "seurat_clusters"
DefaultAssay(PCB) <- "RNA"

new.cluster.ids = c(
  '0' =	'Mature B',
  '1' =	'Mature B',
  '2' =	'Plasma cell',
  '3' =	'Proliferating T',
  '4' =	'Uncertain',
  '5' =	'Proliferating T',
  '6' =	'Plasma cell',
  '7' =	'Uncertain',
  '8' =	'Mature B',
  '9' =	'Doublets',
  '10' =	'Doublets',
  '11' =	'Mature B',
  '12' =	'Doublets',
  '13' =	'Doublets'
)

# Rename the clusters
names(x = new.cluster.ids) <- levels(x = PCB)
PCB <- RenameIdents(object = PCB, new.cluster.ids)

PCB[["celltype_fine"]] <- PCB@active.ident


saveRDS(PCB,"PCBPro1+5_labeled.rds")



print("Script finished!")
