library(Seurat)
library(ggplot2)
library(future)
library(dplyr)
set.seed(10097)
dpi = 500


immune.combined <- readRDS("main_fine_labeled.rds")

# Visualize full UMAP labeled by cell type
png(file="harmony_umap_celltype_fine_rAmbig.png",width=dpi*7,height=dpi*6,units="px",res=dpi)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "celltype_fine",label=TRUE,repel=TRUE,raster=FALSE)
dev.off()

Idents(object=immune.combined) <- "study"

# Visualize UMAP according to severity group
Idents(object=immune.combined) <- "cohort"
png(file="harmony_umap_group_split.png",width=dpi*15,height=dpi*6,units="px",res=dpi)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, split.by = "group",label=FALSE,raster=FALSE)
dev.off()

# Visualize UMAP according to cohort
png(file="harmony_umap_cohort.png",width=dpi*7,height=dpi*6,units="px",res=dpi)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "cohort",label=FALSE,raster=FALSE)
dev.off()

# Visualize UMAP according to cohort
png(file="harmony_umap_cohort_split.png",width=dpi*15,height=dpi*12,units="px",res=dpi)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, split.by = "cohort",label=FALSE,raster=FALSE,ncol=3)
dev.off()

# Visualize UMAP according to BALF vs. PBMC
png(file="harmony_umap_compartment.png",width=dpi*6,height=dpi*6,units="px",res=dpi)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "compartment",label=FALSE,raster=FALSE)
dev.off()

png(file="harmony_umap_compartment_split.png",width=dpi*12,height=dpi*6,units="px",res=dpi)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, split.by = "compartment",label=FALSE,raster=FALSE)
dev.off()

Idents(object=immune.combined) <- "cohort"

liao <- subset(immune.combined, idents = "Liao")
wauters <- subset(immune.combined, idents = "Wauters")
xcontrol <- subset(immune.combined, idents = "BALFxControl")
lee <- subset(immune.combined, idents = "Lee")
arun <- subset(immune.combined, idents = "Arun")
schulte <- subset(immune.combined, idents = "Schulte")
wilk <- subset(immune.combined, idents = "Wilk")

# qc liao
p1 <- VlnPlot(object = liao, features = c("nCount_RNA"),group.by="group",pt.size=0,y.max=80000)
p1
pdf(file="qc_liao_counts_by_group.pdf")
print(p1)
dev.off()

# qc wauters
p1 <- VlnPlot(object = wauters, features = c("nCount_RNA"),group.by="group",pt.size=0,y.max=80000)
p1
pdf(file="qc_wauters_counts_by_group.pdf")
print(p1)
dev.off()

# qc cross control
p1 <- VlnPlot(object = xcontrol, features = c("nCount_RNA"),group.by="group",pt.size=0,y.max=80000)
p1
pdf(file="qc_xcontrol_counts_by_group.pdf")
print(p1)
dev.off()

# qc lee
p1 <- VlnPlot(object = lee, features = c("nCount_RNA"),group.by="group",pt.size=0,y.max=80000)
p1
pdf(file="qc_lee_counts_by_group.pdf")
print(p1)
dev.off()

# qc arun
p1 <- VlnPlot(object = arun, features = c("nCount_RNA"),group.by="group",pt.size=0,y.max=80000)
p1
pdf(file="qc_arun_counts_by_group.pdf")
print(p1)
dev.off()

# qc schulte
p1 <- VlnPlot(object = schulte, features = c("nCount_RNA"),group.by="group",pt.size=0,y.max=80000)
p1
pdf(file="qc_schulte_counts_by_group.pdf")
print(p1)
dev.off()

# qc wilk 
p1 <- VlnPlot(object = wilk, features = c("nCount_RNA"),group.by="group",pt.size=0,y.max=80000)
p1
pdf(file="qc_wilk_counts_by_group.pdf")
print(p1)
dev.off()


# Information to make violins for Lee
# Idents(object=lee) <- "donor"
# print(levels(x=lee))
# levels(x = lee) <- c("Lee-HC1","Lee-HC2","Lee-HC3","Lee-HC4","Lee-M1","Lee-M2","Lee-M3","Lee-S1","Lee-S2","Lee-S3","Lee-S4")
# print(levels(x=lee))
# lee <- StashIdent(object = lee, save.name = 'donor')
# full_obj <- lee
# data_type <- "Lee"

full_obj <- liao
data_type <- "Liao"
Idents(object=liao) <- "donor"
liao_levels <-  c("Liao-HC1","Liao-HC2","Liao-HC3","Liao-M1","Liao-M2","Liao-M3","Liao-S1","Liao-S2","Liao-S3","Liao-S4","Liao-S","Liao-S6")
liao@active.ident <- factor(x = liao@active.ident, levels = liao_levels)
liao[["donor"]] <- liao@active.ident

celltypes <-c("CD14 Mono","cDC2","CD14 Mono", "Inflammatory MP")
degs <- c("CD74","HLA-DPB1","CLU","CLU")

if (data_type=="Wauters"){
  colors=c("#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE
           "#D0E6A5", #86E3CE
           "#D0E6A5", #86E3CE
           "#D0E6A5", #86E3CE
           "#D0E6A5", #86E3CE
           "#D0E6A5", #86E3CE
           "#D0E6A5", #86E3CE
           "#D0E6A5", #86E3CE
           "#D0E6A5", #86E3CE
           "#FFDD94",
           "#FFDD94",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B")
}

if (data_type=="Liao"){
  colors=c("#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#FFDD94",
           "#FFDD94",
           "#FFDD94",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B")
}

if (data_type=="Lee"){
  # Lee-HC1 Lee-HC2 Lee-HC3 Lee-HC4 Lee-ASX1 Lee-M1 Lee-M2 Lee-M3 Lee-S1 Lee-S2 Lee-S2-R Lee-S3 Lee-S3-R Lee-S4 Lee-S4-R
  colors=c("#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE
           "#FFDD94",
           "#FFDD94",
           "#FFDD94",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B")
}

if (data_type=="Arun"){
  # Arun-HC1 Arun-HC2 Arun-HC3 Arun-HC4 Arun-HC5 Arun-M1 Arun-M2 Arun-M3 Arun-S1 Arun-S2 Arun-S3 Arun-S4
  colors=c("#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE
           "#FFDD94",
           "#FFDD94",
           "#FFDD94",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B",
           "#FA897B")
}

if (data_type=="Wilk"){
  colors=c("#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE 
           "#D0E6A5", #86E3CE
           "#D0E6A5", #86E3CE
           "firebrick",
           "#FA897B",
           "#FA897B",
           "firebrick",
           "firebrick",
           "#FA897B",
           "firebrick",
           "#FA897B")
}


Idents(object = full_obj) <- "group"
seurat_obj <- subset(full_obj, idents = c("Severe","Mild","Healthy Control"))


Idents(object=seurat_obj) <- "celltype_fine"
DefaultAssay(seurat_obj) <- "RNA"

i=1
for (cell in celltypes){
  gene=degs[i]
  
  de_subset <- subset(seurat_obj, idents = cell)
  de_subset <- SetIdent(de_subset, value = "donor")
  p0 <- VlnPlot(de_subset, group.by = "donor",pt.size =0.5,features=gene,cols=colors)
  pdf(file = paste0("VIOLIN_",cell,"_",data_type,"_",gene,".pdf"))
  plot(p0)
  dev.off()
  
  i=i+1
}


print("Script finished!")