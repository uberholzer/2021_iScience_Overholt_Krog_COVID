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
# ########################## NKT CLUSTERING ###########################
# # Include 'NKT' and 'Proliferating' coarse celltypes
# 
# Idents(object=immune.combined) <- "celltype_coarse"
# NKT = subset(immune.combined,idents = c("NKT")) # Get NKT clusters
# DefaultAssay(NKT) <- "RNA"
# 
# # Conduct test to choose the number of PCs to include
# pca.test<- NKT %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))
# 
# var.genes <- pca.test@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', pca.test@assays$RNA@var.features)]
# pca.test <- RunPCA(pca.test, features = var.genes, npcs = 100, verbose = FALSE) # Start by calculating 100 PCs
# 
# # Choose number of PCs
# pca.test <- ProjectDim(object = pca.test)
# pdf(file="NKT1+5_elbowplot.pdf")
# ElbowPlot(object = pca.test,ndims = 100)
# dev.off()
# 
# pct <- pca.test[["pca"]]@stdev / sum(pca.test[["pca"]]@stdev) * 100
# nPCs <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > PCTthresh), decreasing = T)[1] + 1 + n_beyond
# 
# cat("/nNumber of PCs used for NKT: ")
# print(nPCs)
# 
# # Run the actual PCA using the number of PCs calculated
# NKT <- NKT %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))
# 
# removed.genes <- NKT@assays$RNA@var.features[grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', NKT@assays$RNA@var.features)]
# print(paste("Removed genes: ",removed.genes))
# 
# var.genes <- NKT@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', NKT@assays$RNA@var.features)]
# NKT <- RunPCA(NKT, features = var.genes, npcs = nPCs, verbose = FALSE)
# 
# print("Scale and PCA successful!")
# 
# # Run Harmony
# NKT <- NKT %>% RunHarmony("orig.ident", plot_convergence = FALSE) # make sure to specify what variable(s) to regress out here
# 
# # Clustering
# NKT <- NKT %>%
#   RunUMAP(reduction = "harmony", dims = 1:nPCs) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:nPCs) %>%
#   FindClusters(resolution = 0.5) %>%
#   identity()
# 
# print("\nTotal number of cells: ")
# print(ncol(NKT@assays[["RNA"]]))
# 
# # Visualization
# pdf(file="NKT1+5_clusters.pdf")
# DimPlot(NKT, reduction = "umap", label = TRUE, pt.size = .1)
# dev.off()
# 
# pdf(file="NKT1+5_donors.pdf",width=18,height=6)
# DimPlot(object = NKT, reduction = "umap", pt.size = .1, group.by = "orig.ident",label=FALSE)
# dev.off()
# 
# pdf(file="NKT1+5_study.pdf",width=6,height=6)
# DimPlot(object = NKT, reduction = "umap", pt.size = .1, group.by = "study",label=FALSE)
# dev.off()
# 
# saveRDS(NKT,"NKT1+5.rds")



########################## MARKER DOTPLOTS AND HEATMAPS ###########################


NKT <- readRDS("NKT1+5.rds")

# Idents(object=NKT) <- "seurat_clusters"
# DefaultAssay(NKT) <- "RNA"

# # Make heatmap
# # First scale all genes to make a complete heatmap
# NKT <- ScaleData(NKT, features = rownames(NKT))
#
# # # Parallelize FindMarkers step
# plan("multiprocess") # works with 4 workers
# plan()
# options(future.globals.maxSize = 10000 * 1024^9,seed=TRUE)
#
# # find markers for every cluster compared to all remaining cells, report only the positive ones
# NKT@misc$markers <- FindAllMarkers(object = NKT, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
# write.table(NKT@misc$markers,file='NKT-markers_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')
#
# ##add average expression information
# # NKT@misc$averageExpression = AverageExpression(object = NKT)
# # write.table(NKT@misc$averageExpression$RNA,file='NKT-average_MAST.txt',row.names = TRUE,quote = FALSE,sep = '\t')
#
# hc.markers = read.delim2("NKT1+5rPro05_markers.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
# hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
# 
# tt = DoHeatmap(object = NKT, features = top10$gene) + NoLegend()  +
#   scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdPu"))
# ggplot2::ggsave(file="NKT1+5rPro05_heatmap.pdf",plot = tt,device = 'pdf',width = 16, height = 12, units = "in",dpi = 300)


Idents(object=NKT) <- "seurat_clusters"
DefaultAssay(NKT) <- "RNA"

# # Dotplots for cell type annotation
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

p4 <- DotPlot(NKT, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "NKT_dotplot_simple.pdf",height=10,width=15)
plot(p4)
p4
dev.off()

features=rev(c("EMP2","TPPP3","KRT18","KRT15","SCGB1A1","TPSB2","TPSAB1","TYMS","MKI67","IGHA1","IGHG4","IGHM",
               "IGHD","CD22","CD19","MS4A1","CD79A","TRAV1-2","CXCR6","GZMK","KLRB1","EOMES","KLRF1","GZMB","GNLY","KLRD1","KLRG1","LIM2","CD28",
               "CD69","FAS","FOXP3","CTLA4","TRDV3","TRDV2","TRDV1","TRGC1","CXCR3","IL2RA","IL2RB","CD58","ITGAL","CD44","IL7R","IL2RG","TCF7","CCR7","SELL","CD8A","CD4","CD3D",
               "CXCR2","FCGR3B","HCAR3","PI3","ELANE","LILRA4","IL3RA","THBD","CD1C","CLEC9A","XCR1","FABP4","SPP1",
               "FCGR3A","CD14","S100A8","FCN1","CD68","PPBP","HBB"))

p4 <- DotPlot(NKT, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "NKT1+5_dotplot.pdf",height=10,width=42)
plot(p4)
p4
dev.off()


########################## CLUSTER ANNOTATION ###########################

Idents(object=NKT) <- "seurat_clusters"
DefaultAssay(NKT) <- "RNA"

new.cluster.ids = c(
  '0' = 	'NK',
  '1' = 	'CD4 Naive T',
  '2' = 	'CD8 Effector T',
  '3' = 	'CD4 Naive T',
  '4' = 	'CD8 Effector T',
  '5' = 	'CD4 Treg',
  '6' = 	'CD8 Naive T',
  '7' = 	'CD8 Effector T',
  '8' = 	'Doublets',
  '9' = 	'Doublets',
  '10' = 	'CD8 Effector T',
  '11' = 	'CD4 Naive T',
  '12' = 	'Uncertain',
  '13' = 	'Uncertain'
)

# Rename the clusters
names(x = new.cluster.ids) <- levels(x = NKT)
NKT <- RenameIdents(object = NKT, new.cluster.ids)

NKT[["celltype_fine"]] <- NKT@active.ident


Idents(object = NKT) <- "seurat_clusters"

NKT[['percent.igtot']] <- (PercentageFeatureSet(NKT, pattern = "^IGK", assay = 'RNA')
                                +PercentageFeatureSet(NKT, pattern = "^IGH", assay = 'RNA')
                                +PercentageFeatureSet(NKT, pattern = "^IGL", assay = 'RNA'))
# Plot original
p0 <- VlnPlot(NKT,pt.size =0.5,features="percent.igtot")
pdf(file = "NKT_pctIgtot.pdf")
plot(p0)
dev.off()

# Filter cells by total percentage of immunoglobulin and BCR genes
NKT <- subset(x = NKT, subset = percent.igtot < 5)

# Plot original
p0 <- VlnPlot(NKT,pt.size =0.5,features="percent.igtot")
pdf(file = "NKT_pctIgtot_filt.pdf")
plot(p0)
dev.off()


saveRDS(NKT,"NKT1+5_labeled.rds")

print("Script finished!")
