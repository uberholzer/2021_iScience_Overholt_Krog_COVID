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
# ########################## DENDRITIC CLUSTERING ###########################
# # Include 'Dendritic' coarse annotations only
# 
# Idents(object=immune.combined) <- "celltype_coarse"
# dendritic = subset(immune.combined,idents = c("Dendritic")) # Get Dendritic cell clusters
# DefaultAssay(dendritic) <- "RNA"
# 
# # Conduct test to choose the number of PCs to include
# pca.test<- dendritic %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))
# 
# var.genes <- pca.test@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', pca.test@assays$RNA@var.features)]
# pca.test <- RunPCA(pca.test, features = var.genes, npcs = 100, verbose = FALSE) # Start by calculating 100 PCs
# 
# # Choose number of PCs
# pca.test <- ProjectDim(object = pca.test)
# pdf(file="dendritic1+5_elbowplot.pdf")
# ElbowPlot(object = pca.test,ndims = 100)
# dev.off()
# 
# pct <- pca.test[["pca"]]@stdev / sum(pca.test[["pca"]]@stdev) * 100
# nPCs <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > PCTthresh), decreasing = T)[1] + 1 + n_beyond
# 
# cat("/nNumber of PCs used for dendritic: ")
# print(nPCs)
# 
# # Run the actual PCA using the number of PCs calculated
# dendritic <- dendritic %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))
# 
# removed.genes <- dendritic@assays$RNA@var.features[grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', dendritic@assays$RNA@var.features)]
# print(paste("Removed genes: ",removed.genes))
# 
# var.genes <- dendritic@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', dendritic@assays$RNA@var.features)]
# dendritic <- RunPCA(dendritic, features = var.genes, npcs = nPCs, verbose = FALSE)
# 
# print("Scale and PCA successful!")
# 
# # Run Harmony
# dendritic <- dendritic %>% RunHarmony("orig.ident", plot_convergence = FALSE) # make sure to specify what variable(s) to regress out here
# 
# # # Clustering
# dendritic <- dendritic %>%
#   RunUMAP(reduction = "harmony", dims = 1:nPCs) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:nPCs) %>%
#   FindClusters(resolution = 0.5) %>%
#   identity()
# 
# print("\nTotal number of cells: ")
# print(ncol(dendritic@assays[["RNA"]]))
# 
# # Visualization
# pdf(file="dendritic1+5_clusters.pdf")
# DimPlot(dendritic, reduction = "umap", label = TRUE, pt.size = .1)
# dev.off()
# 
# pdf(file="dendritic1+5_donors.pdf",width=18,height=6)
# DimPlot(object = dendritic, reduction = "umap", pt.size = .1, group.by = "orig.ident",label=FALSE)
# dev.off()
# 
# pdf(file="dendritic1+5_study.pdf",width=6,height=6)
# DimPlot(object = dendritic, reduction = "umap", pt.size = .1, group.by = "study",label=FALSE)
# dev.off()
# 
# saveRDS(dendritic,"dendritic1+5.rds")

# # find markers for every cluster compared to all remaining cells, report only the positive ones
# dendritic@misc$markers <- FindAllMarkers(object = dendritic, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
# write.table(dendritic@misc$markers,file='dendritic-markers_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')
# ##add average expression information
#
# dendritic@misc$averageExpression = AverageExpression(object = dendritic)
# write.table(dendritic@misc$averageExpression$RNA,file='dendritic-average_MAST.txt',row.names = TRUE,quote = FALSE,sep = '\t')
#
#
# Idents(object=dendritic) <- "seurat_clusters"
#
# hc.markers = read.delim2("dendritic-markers_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
# hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
#
# tt = DoHeatmap(object = dendritic, features = top10$gene) + NoLegend()  +
#   scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdPu"))
# ggplot2::ggsave(file="dendritic-heatmap.pdf",plot = tt,device = 'pdf',width = 16, height = 12, units = "in",dpi = dpi)



########################## MARKER DOTPLOTS AND HEATMAPS ###########################
dendritic <- readRDS("dendritic1+5.rds")

Idents(object=dendritic) <- "seurat_clusters"
DefaultAssay(dendritic) <- "RNA"

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

p4 <- DotPlot(dendritic, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "dendritic_dotplot_simple.pdf",height=10,width=15)
plot(p4)
p4
dev.off()

features=rev(c("EMP2","TPPP3","KRT18","KRT15","SCGB1A1","TPSB2","TPSAB1","TYMS","MKI67","IGHA1","IGHG4","IGHM",
               "IGHD","CD22","CD19","MS4A1","CD79A","CXCR6","GZMK","KLRB1","EOMES","KLRF1","GZMB","GNLY","KLRD1",
               "CD69","FAS","FOXP3","CTLA4","IL2RA","CD44","IL7R","IL2RG","TCF7","CCR7","SELL","CD8A","CD4","CD3D",
               "CXCR2","FCGR3B","HCAR3","PI3","ELANE", "LILRA4","IL3RA","THBD","AXL","SIRPA","IRF4","CD1C","CLEC9A","XCR1","FABP4","SPP1",
               "FCGR3A","CD14","S100A8","FCN1","CD68","PPBP","HBB"))

p4 <- DotPlot(dendritic, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "dendritic1+5_dotplot.pdf",height=10,width=30)
plot(p4)
p4
dev.off()

########################## CLUSTER ANNOTATION ###########################

Idents(object=dendritic) <- "seurat_clusters"
DefaultAssay(dendritic) <- "RNA"

new.cluster.ids = c(
  '0' =	'pDC',
  '1' =	'cDC2',
  '2' =	'cDC2',
  '3' =	'Uncertain',
  '4' =	'Doublets',
  '5' =	'Uncertain',
  '6' =	'Uncertain',
  '7' =	'Uncertain',
  '8' =	'Doublets',
  '9' =	'cDC2',
  '10' =	'cDC2',
  '11' =	'cDC1',
  '12' =	'Uncertain',
  '13' =	'Doublets',
  '14' =	'pDC'
)

# Rename the clusters
names(x = new.cluster.ids) <- levels(x = dendritic)
dendritic <- RenameIdents(object = dendritic, new.cluster.ids)

dendritic[["celltype_fine"]] <- dendritic@active.ident


Idents(object = dendritic) <- "seurat_clusters"

dendritic[['percent.igtot']] <- (PercentageFeatureSet(dendritic, pattern = "^IGK", assay = 'RNA')
                               +PercentageFeatureSet(dendritic, pattern = "^IGH", assay = 'RNA')
                               +PercentageFeatureSet(dendritic, pattern = "^IGL", assay = 'RNA'))
# Plot original
p0 <- VlnPlot(dendritic,pt.size =0.5,features="percent.igtot")
pdf(file = "dendritic_pctIgtot.pdf")
plot(p0)
dev.off()

# Filter cells by total percentage of immunoglobulin and BCR genes
dendritic <- subset(x = dendritic, subset = percent.igtot < 5)

# Plot original
p0 <- VlnPlot(dendritic,pt.size =0.5,features="percent.igtot")
pdf(file = "dendritic_pctIgtot_filt.pdf")
plot(p0)
dev.off()


saveRDS(dendritic,"dendritic1+5_labeled.rds")

print("Script finished!")
