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
# n_beyond <- 0
# 
# 
# ########################## MYELOID CLUSTERING ###########################
# # Include 'Myeloid' coarse celltypes only
# 
# Idents(object=immune.combined) <- "celltype_coarse"
# myeloid = subset(immune.combined,idents = c("Myeloid")) # Get the Myeloid clusters
# DefaultAssay(myeloid) <- "RNA"
# 
# # Conduct test to choose the number of PCs to include
# pca.test<- myeloid %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))
# 
# var.genes <- pca.test@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', pca.test@assays$RNA@var.features)]
# pca.test <- RunPCA(pca.test, features = var.genes, npcs = 100, verbose = FALSE) # Start by calculating 100 PCs
# 
# # Choose number of PCs
# pca.test <- ProjectDim(object = pca.test)
# pdf(file="myeloid1+0_elbowplot.pdf")
# ElbowPlot(object = pca.test,ndims = 100)
# dev.off()
# 
# pct <- pca.test[["pca"]]@stdev / sum(pca.test[["pca"]]@stdev) * 100
# nPCs <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > PCTthresh), decreasing = T)[1] + 1 + n_beyond
# 
# cat("/nNumber of PCs used for myeloid: ")
# print(nPCs)
# 
# # Run the actual PCA using the number of PCs calculated
# myeloid <- myeloid %>%
#   Seurat::NormalizeData(verbose = FALSE) %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))
# 
# removed.genes <- myeloid@assays$RNA@var.features[grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', myeloid@assays$RNA@var.features)]
# print(paste("Removed genes: ",removed.genes))
# 
# var.genes <- myeloid@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', myeloid@assays$RNA@var.features)]
# myeloid <- RunPCA(myeloid, features = var.genes, npcs = nPCs, verbose = FALSE)
# 
# print("Scale and PCA successful!")
# 
# # Run Harmony
# myeloid <- myeloid %>% RunHarmony("orig.ident", plot_convergence = FALSE) # make sure to specify what variable(s) to regress out here
# 
# # Clustering
# myeloid <- myeloid %>%
#   RunUMAP(reduction = "harmony", dims = 1:nPCs) %>%
#   FindNeighbors(reduction = "harmony", dims = 1:nPCs) %>%
#   FindClusters(resolution = 0.5) %>%
#   identity()
# 
# print("\nTotal number of cells: ")
# print(ncol(myeloid@assays[["RNA"]]))
# 
# # Visualization
# pdf(file="myeloid1+0_clusters.pdf")
# DimPlot(myeloid, reduction = "umap", label = TRUE, pt.size = .1)
# dev.off()
# 
# pdf(file="myeloid1+0_donors.pdf",width=18,height=6)
# DimPlot(object = myeloid, reduction = "umap", pt.size = .1, group.by = "orig.ident",label=FALSE)
# dev.off()
# 
# pdf(file="myeloid1+0_study.pdf",width=6,height=6)
# DimPlot(object = myeloid, reduction = "umap", pt.size = .1, group.by = "study",label=FALSE)
# dev.off()
# 
# saveRDS(myeloid,"myeloid1+0.rds")


# ########################## MARKER DOTPLOTS AND HEATMAPS ###########################

myeloid <- readRDS("myeloid1+0.rds")

Idents(object=myeloid) <- "seurat_clusters"
DefaultAssay(myeloid) <- "RNA"
#
# # # Dotplots for cell type annotation
# # Dotplots for cell type annotation
features=rev(c("TPPP3","SCGB1A1", # Epithelial
               "TYMS","MKI67", # Proliferating
               "IGHA1","IGHG4", # Plasma cells
               "IGHM","IGHD","MS4A1", # B cells
               "NKG7","KLRD1", # NK/effector
               "CD8A","CD4","CD3D", # T cells
               "LILRA4","IL3RA","CD1C", # Dendritic cells
               "IL1B","FCGR3A","CD14","CD68","CD16","VCAN", # Monocyte/macrophages
               "PPBP", # Platelets
               "HBB")) # RBCs

p4 <- DotPlot(myeloid, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "myeloid1+0_dotplot_simple.pdf",height=10,width=15)
plot(p4)
p4
dev.off()

features=rev(c("EMP2","TPPP3","KRT18","KRT15","SCGB1A1","TPSB2","TPSAB1","TYMS","MKI67","IGHA1","IGHG4","IGHM",
               "IGHD","CD22","CD19","MS4A1","CD79A","CXCR6","GZMK","KLRB1","EOMES","KLRF1","GZMB","GNLY","KLRD1",
               "CD69","FAS","FOXP3","CTLA4","IL2RA","CD44","IL7R","IL2RG","TCF7","CCR7","SELL","CD8A","CD4","CD3D",
               "CEACAM8","CXCR2","FCGR3B","HCAR3","PI3","ELANE", "LILRA4","IL3RA","THBD","CD1C","CLEC9A","XCR1","IL1B","IL1RN","TNF","CXCL2","CCL3","CXCL3","CCL4","FABP4","SPP1",
               "FCGR3A","CD14","S100A8","FCN1","CD68","ITGAM","MSR1","MRC1","SIGLEC1","PPBP","HBB"))

p4 <- DotPlot(myeloid, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "myeloid1+0_dotplot.pdf",height=10,width=40)
plot(p4)
p4
dev.off()

# Make heatmap
# First scale all genes to make a complete heatmap
# myeloid <- ScaleData(myeloid, features = rownames(myeloid))

# # Parallelize FindMarkers step
# plan("multiprocess") # works with 4 workers
# plan()
# options(future.globals.maxSize = 10000 * 1024^9,seed=TRUE)
# 
# # find markers for every cluster compared to all remaining cells, report only the positive ones
# myeloid@misc$markers <- FindAllMarkers(object = myeloid, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
# write.table(myeloid@misc$markers,file='myeloid1+0_markers.txt',row.names = FALSE,quote = FALSE,sep = '\t')

##add average expression information
# myeloid@misc$averageExpression = AverageExpression(object = myeloid)
# write.table(myeloid@misc$averageExpression$RNA,file='myeloid-average_MAST.txt',row.names = TRUE,quote = FALSE,sep = '\t')

# hc.markers = read.delim2("myeloid1+0_markers.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
# hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
# 
# tt = DoHeatmap(object = myeloid, features = top10$gene) + NoLegend()  +
#   scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdPu"))
# ggplot2::ggsave(file="myeloid1+0_heatmap.pdf",plot = tt,device = 'pdf',width = 16, height = 12, units = "in",dpi = 300)

########################## CLUSTER ANNOTATION ###########################

Idents(object=myeloid) <- "seurat_clusters"
DefaultAssay(myeloid) <- "RNA"

new.cluster.ids = c(
  '0' =	  'AM',
  '1' =	  'CD14 Mono',
  '2' =	  'AM',
  '3' =	  'AM',
  '4' =	  'Inflammatory MP',
  '5' =	  'AM',
  '6' =	  'AM',
  '7' =	  'CD16 Mono',
  '8' =	  'Neutrophil',
  '9' =	  'AM',
  '10' =	'Doublets',
  '11' =	'AM',
  '12' =	'Doublets',
  '13' =	'Intermediate Mono',
  '14' =	'Doublets',
  '15' =	'Doublets',
  '16' =	'Doublets',
  '17' =	'Doublets'
)

# Rename the clusters
names(x = new.cluster.ids) <- levels(x = myeloid)
myeloid <- RenameIdents(object = myeloid, new.cluster.ids)

myeloid[["celltype_fine"]] <- myeloid@active.ident



Idents(object = myeloid) <- "seurat_clusters"

myeloid[['percent.igtot']] <- (PercentageFeatureSet(myeloid, pattern = "^IGK", assay = 'RNA')
                           +PercentageFeatureSet(myeloid, pattern = "^IGH", assay = 'RNA')
                           +PercentageFeatureSet(myeloid, pattern = "^IGL", assay = 'RNA'))
# Plot original
p0 <- VlnPlot(myeloid,pt.size =0.5,features="percent.igtot")
pdf(file = "myeloid_pctIgtot.pdf")
plot(p0)
dev.off()

# Filter cells by total percentage of immunoglobulin and BCR genes
myeloid <- subset(x = myeloid, subset = percent.igtot < 5)

# Plot original
p0 <- VlnPlot(myeloid,pt.size =0.5,features="percent.igtot")
pdf(file = "myeloid_pctIgtot_filt.pdf")
plot(p0)
dev.off()


saveRDS(myeloid,"myeloid1+0_labeled.rds")


print("Script finished!")
