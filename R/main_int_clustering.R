library(Seurat)
library(ggplot2)
library(future)
library(dplyr)
library(harmony)
set.seed(10097)

# CLUSTER ANNOTATION AND SUBCLUSTERING 

immune.combined <- readRDS("main_harmony_integrated.rds")

# Set clustering paramters
PCTthresh <- 0.05 # Set threshold here for determining the number of PCs to use
n_beyond <- 5

########################## MARKER DOTPLOTS AND HEATMAPS ###########################

Idents(object=immune.combined) <- "seurat_clusters"
DefaultAssay(immune.combined) <- "RNA" 

# Dotplots for cell type annotation
features=rev(c("TPPP3","KRT18","SCGB1A1", # Epithelial
               "TYMS","MKI67", # Proliferating
               "IGHA1","IGHG4", # Plasma cells
               "IGHM","IGHD","MS4A1", # B cells
               "NKG7","KLRD1", # NK/effector
               "CD8A","CD4","CD3D", # T cells
               "LILRA4","IL3RA","CD1C", # Dendritic cells
               "FCGR3A","CD14","CD68","CD16","VCAN", # Monocyte/macrophages
               "PPBP", # Platelets
               "SNCA","HBB")) # RBCs

p4 <- DotPlot(immune.combined, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "main_dotplot.pdf",height=10,width=15)
plot(p4)
p4
dev.off()

features=rev(c("EMP2","TPPP3","KRT18","KRT15","SCGB1A1","TPSB2","TPSAB1","TYMS","MKI67","IGHA1","IGHG4","IGHM",
               "IGHD","CD22","CD19","MS4A1","CD79A","CXCR6","GZMK","KLRB1","EOMES","KLRF1","GZMB","GNLY","KLRD1",
               "CD69","FAS","FOXP3","CTLA4","IL2RA","CD44","IL7R","IL2RG","TCF7","CCR7","SELL","CD8A","CD4","CD3D",
               "CXCR2","FCGR3B","HCAR3","PI3","ELANE", "LILRA4","IL3RA","THBD","CD1C","CLEC9A","XCR1","FABP4","SPP1",
               "FCGR3A","CD14","S100A8","FCN1","CD68","PPBP")) 

p4 <- DotPlot(immune.combined, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "main_dotplot_full.pdf",height=10,width=30)
plot(p4)
p4
dev.off()

# Calculate markers for the Seurat clusters using MAST (computationally expensive step)
# find markers for every cluster compared to all remaining cells, report only the positive ones
# Parallelize

# plan("multiprocess") # works with 4 workers
# plan()
# options(future.globals.maxSize = 10000 * 1024^9,seed=TRUE) #TRUE here?
# 
# immune.combined@misc$markers <- FindAllMarkers(object = immune.combined, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
# print("Markers found!")
# write.table(immune.combined@misc$markers,file='main_marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')
# 
# hc.markers = read.delim2("main_marker_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
# hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
# 
# write.csv(top10,"main_top10_markers.csv", row.names = FALSE) # Write the markers to a file
# 
# # Draw heatmap of marker genes
# tt1 = DoHeatmap(object = subset(immune.combined, downsample = 500), features = top10$gene)  + 
#   scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdPu"))
# # + NoLegend() # if needed
# ggplot2::ggsave(file="main_marker_heatmap.pdf",plot = tt1,device = 'pdf',width = 24, height = 16, units = "in",dpi = 300,limitsize = FALSE)
# 

########################## COARSE CLUSTERING ###########################

new.cluster.ids = c(
  
  '0'='NKT',
  '1'='Myeloid',
  '2'='Myeloid',
  '3'='NKT',
  '4'='Myeloid',
  '5'='B cell',
  '6'='Dendritic',
  '7'='Proliferating',
  '8'='Myeloid',
  '9'='Plasma cell',
  '10'='Platelet',
  '11'='Epithelial',
  '12'='Epithelial',
  '13'='Myeloid',
  '14'='Dendritic',
  '15'='Erythrocyte',
  '16'='Myeloid',
  '17'='NKT',
  '18'='Myeloid',
  '19'='Myeloid',
  '20'='NKT',
  '21'='NKT'
  
)

# Rename the clusters
names(x = new.cluster.ids) <- levels(x = immune.combined)
immune.combined <- RenameIdents(object = immune.combined, new.cluster.ids)

immune.combined[["celltype_coarse"]] <- immune.combined@active.ident

saveRDS(immune.combined,"main_coarse_labeled.rds")


print("Script finished!")
