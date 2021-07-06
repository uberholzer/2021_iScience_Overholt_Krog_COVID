library(Seurat)
library(ggplot2)
library(future)
library(dplyr)
set.seed(10097)

# Read in coarse annotations
immune.combined <- readRDS("main_coarse_labeled.rds")

########################### TRANSFER METADATA FROM SUBSETTED OBJECTS ###########################
# Load subsets
dendriticObj <- readRDS("dendritic1+5_labeled.rds")
myeloidObj <- readRDS("myeloid1+0_labeled.rds")
NKTObj <- readRDS("NKT1+5_labeled.rds")
PCBObj <- readRDS("PCBPro1+5_labeled.rds")

# Create metadata dataframe to edit
fine_annot <- immune.combined@meta.data['celltype_coarse']
names(fine_annot)[names(fine_annot) == "celltype_coarse"] <- "celltype_fine"

# Unfactorize
fine_annot$celltype_fine <- as.character(fine_annot$celltype_fine)

# Take out the coarse annotations we are replacing
fine_annot$celltype_fine[fine_annot$celltype_fine=="Myeloid"] <- "UNASSIGNED"
fine_annot$celltype_fine[fine_annot$celltype_fine=="Dendritic"] <- "UNASSIGNED"
fine_annot$celltype_fine[fine_annot$celltype_fine=="NKT"] <- "UNASSIGNED"
fine_annot$celltype_fine[fine_annot$celltype_fine=="B cell"] <- "UNASSIGNED"
fine_annot$celltype_fine[fine_annot$celltype_fine=="Plasma cell"] <- "UNASSIGNED"
fine_annot$celltype_fine[fine_annot$celltype_fine=="Proliferating"] <- "UNASSIGNED"

# Add dendritic annotations
dendritic <- dendriticObj@meta.data['celltype_fine']
dendritic$celltype_fine <- as.character(dendritic$celltype_fine)
fine_annot[match(row.names(dendritic), row.names(fine_annot)), ] <- dendritic # Replace annotations in fine_annot with annotations in dendritic for all matching barcodes
# This method was shown to work in a toy script called 'recombine_tester.Rmd' that was saved locally

# Add PCB annotations
PCB <- PCBObj@meta.data['celltype_fine']
PCB$celltype_fine <- as.character(PCB$celltype_fine)
fine_annot[match(row.names(PCB), row.names(fine_annot)), ] <- PCB # Replace annotations in fine_annot with annotations in PCB for all matching barcodes

# Add myeloid annotations
myeloid <- myeloidObj@meta.data['celltype_fine']
myeloid$celltype_fine <- as.character(myeloid$celltype_fine)
fine_annot[match(row.names(myeloid), row.names(fine_annot)), ] <- myeloid # Replace annotations in fine_annot with annotations in myeloid for all matching barcodes

# Add NKT annotations
NKT <- NKTObj@meta.data['celltype_fine']
NKT$celltype_fine <- as.character(NKT$celltype_fine)
fine_annot[match(row.names(NKT), row.names(fine_annot)), ] <- NKT # Replace annotations in fine_annot with annotations in NKT for all matching barcodes

# Re-factorize the full fine annotation dataframe
fine_annot$celltype_fine <- as.factor(fine_annot$celltype_fine)

# Add fine annotation metadata
immune.combined <- AddMetaData(
  object = immune.combined,
  metadata = fine_annot,
  col.name = 'celltype_fine'
)

# DELETE cells that are still unassigned as a result of being removed during subclustering
immune.combined <- subset(immune.combined, subset = celltype_fine != "UNASSIGNED")

# Visualize
pdf(file="harmony_umap_celltype_coarse.pdf",width=7,height=6)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "celltype_coarse",label=TRUE)
dev.off()

pdf(file="harmony_umap_celltype_fine.pdf",width=9,height=6)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "celltype_fine",label=TRUE,repel=TRUE)
dev.off()

# Eliminate the Doublets and Uncertain clusters
immune.combined <- subset(x=immune.combined, subset = celltype_fine != "Doublets" & celltype_fine != "Uncertain")

# Visualize again
pdf(file="harmony_umap_celltype_fine_rAmbig.pdf",width=7,height=6)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "celltype_fine",label=TRUE,repel=TRUE,raster=FALSE)
dev.off()

# Add donor information
# See "all_donors" table
Idents(object=immune.combined) <- "orig.ident"
immune.combined <- RenameIdents(object = immune.combined,
                                'C51'='Liao-HC1',
                                'C52'='Liao-HC2',
                                'C100'='Liao-HC3',
                                'C141'='Liao-M1',
                                'C142'='Liao-M2',
                                'C144'='Liao-M3',
                                'C145'='Liao-S1',
                                'C143'='Liao-S2',
                                'C146'='Liao-S3',
                                'C148'='Liao-S4',
                                'C149'='Liao-S5',
                                'C152'='Liao-S6',
                                
                                'GSM3660650'='Tabib-HC1',
                                
                                'nCoV 1 scRNA-seq' = 'Lee-S1', # Patient C1
                                'nCoV 2 scRNA-seq' = 'Lee-M1', # Patient C2
                                'nCoV 3 scRNA-seq' = 'Lee-S2',  # Patient C3
                                'nCoV 5 scRNA-seq' = 'Lee-M2', # Patient C4
                                'nCoV 6 scRNA-seq' = 'Lee-M3', # Patient C5
                                'nCoV 7 scRNA-seq' = 'Lee-S3', # Patient C6
                                'nCoV 9 scRNA-seq' = 'Lee-S4', # Patient C7
                                'Normal 1 scRNA-seq' = 'Lee-HC1', 
                                'Normal 2 scRNA-seq' = 'Lee-HC2', 
                                'Normal 3 scRNA-seq' = 'Lee-HC3',
                                'Normal 4 scRNA-seq' = 'Lee-HC4',
                                
                                'cov1'  = 'Arun-S1',
                                'cov2'  = 'Arun-M1', 
                                'cov3'  = 'Arun-M2', 
                                'cov4'  = 'Arun-S2',
                                'cov7'  = 'Arun-HC1',
                                'cov8'  = 'Arun-HC2',
                                'cov9'  = 'Arun-HC3',
                                'cov10' = 'Arun-S3',
                                'cov11' = 'Arun-S4',
                                'cov12' = 'Arun-M3',
                                'cov17' = 'Arun-HC4',
                                'cov18' = 'Arun-HC5',
                                
                                'covid_555_1'='Wilk-S1-N', # Non-ARDS
                                'covid_556'='Wilk-S2-N', # Non-ARDS
                                'covid_557'='Wilk-S3-A', # ARDS
                                'covid_558'='Wilk-S4-A', # ARDS
                                'covid_559'='Wilk-S5-A', # ARDS
                                'covid_560'='Wilk-S6-A', # ARDS
                                'covid_561'='Wilk-S7-N', # Non-ARDS
                                'HIP002'='Wilk-HC1',
                                'HIP015'='Wilk-HC2',
                                'HIP023'='Wilk-HC3',
                                'HIP043'='Wilk-HC4',
                                'HIP044'='Wilk-HC5',
                                'HIP045'='Wilk-HC6',
                                
                                'C19-CB-0001'='Schulte-M1',
                                'C19-CB-0002'='Schulte-M2',
                                'C19-CB-0003'='Schulte-M3',
                                'C19-CB-0204'='Schulte-M4',
                                'C19-CB-0214'='Schulte-M5',
                                'C19-CB-0008'='Schulte-S1',
                                'C19-CB-0009'='Schulte-S2',
                                'C19-CB-0011'='Schulte-S3',
                                'C19-CB-0012'='Schulte-S4',
                                'C19-CB-0013'='Schulte-S5',
                                'C19-CB-0016'='Schulte-S6',
                                'C19-CB-0020'='Schulte-S7',
                                'C19-CB-0021'='Schulte-S8',
                                'C19-CB-0198'='Schulte-S9',
                                'C19-CB-0199'='Schulte-S10',
                                
                                'BAL009' =	'Wauters-M1',
                                'BAL012' =	'Wauters-S1',
                                'BAL013' =	'Wauters-S2',
                                'BAL014' =	'Wauters-S3',
                                'BAL015' =	'Wauters-S4',
                                'BAL016' =	'Wauters-S5',
                                'BAL022' =	'Wauters-S6',
                                'BAL023' =	'Wauters-S7',
                                'BAL024' =	'Wauters-S8',
                                'BAL025' =	'Wauters-S9',
                                'BAL032' =	'Wauters-S10',
                                'BAL035' =	'Wauters-S11',
                                'BAL039' =	'Wauters-S12',
                                'BAL040' =	'Wauters-S13',
                                'BAL027' =	'Wauters-S14',
                                'BAL031' =	'Wauters-S15',
                                'BAL037' =	'Wauters-M2',
                                'BAL021' =	'Wauters-S16',
                                'BAL026' =	'Wauters-S17',
                                'BAL033' =	'Wauters-S18',
                                'BAL034' =	'Wauters-S19',
                                'BAL020' =	'Wauters-S20',
                                
                                'GSM4593897_sample_10_UMI_counts.csv'='Mould-HC10',
                                'GSM4593896_sample_9_UMI_counts.csv'='Mould-HC9',
                                'GSM4593895_sample_8_UMI_counts.csv'='Mould-HC8',
                                'GSM4593894_sample_7_UMI_counts.csv'='Mould-HC7',
                                'GSM4593893_sample_6_UMI_counts.csv'='Mould-HC6',
                                'GSM4593892_sample_5_UMI_counts.csv'='Mould-HC5',
                                'GSM4593891_sample_4_UMI_counts.csv'='Mould-HC4',
                                'GSM4593890_sample_3_UMI_counts.csv'='Mould-HC3',
                                'GSM4593889_sample_2_UMI_counts.csv'='Mould-HC2',
                                'GSM4593888_sample_1_UMI_counts.csv'='Mould-HC1')

immune.combined[["donor"]] <- immune.combined@active.ident

# Add severity group information
Idents(object=immune.combined) <- "orig.ident"
immune.combined <- RenameIdents(object = immune.combined,
                                'C51'='Healthy Control',
                                'C52'='Healthy Control',
                                'C100'='Healthy Control',
                                'C141'='Mild',
                                'C142'='Mild',
                                'C144'='Mild',
                                'C145'='Severe',
                                'C143'='Severe',
                                'C146'='Severe',
                                'C148'='Severe',
                                'C149'='Severe',
                                'C152'='Severe',
                                
                                'GSM3660650'='Healthy Control',
                                
                                'nCoV 1 scRNA-seq'  = 'Severe', # Patient C1
                                'nCoV 2 scRNA-seq'  = 'Mild', # Patient C2
                                'nCoV 3 scRNA-seq'  = 'Severe',  # Patient C3
                                'nCoV 5 scRNA-seq'  = 'Mild', # Patient C4
                                'nCoV 6 scRNA-seq'  = 'Mild', # Patient C5
                                'nCoV 7 scRNA-seq'  = 'Severe', # Patient C6
                                'nCoV 9 scRNA-seq'  = 'Severe', # Patient C7
                                'Normal 1 scRNA-seq' = 'Healthy Control', 
                                'Normal 2 scRNA-seq' = 'Healthy Control', 
                                'Normal 3 scRNA-seq' = 'Healthy Control',
                                'Normal 4 scRNA-seq' = 'Healthy Control',
                                
                                'cov1'  = 'Severe',
                                'cov2'  = 'Mild', 
                                'cov3'  = 'Mild', 
                                'cov4'  = 'Severe',
                                'cov7'  = 'Healthy Control',
                                'cov8'  = 'Healthy Control',
                                'cov9'  = 'Healthy Control',
                                'cov10' = 'Severe',
                                'cov11' = 'Severe',
                                'cov12' = 'Mild',
                                'cov17' = 'Healthy Control',
                                'cov18' = 'Healthy Control',
                                
                                'covid_555_1'='Severe', # Non-ARDS
                                'covid_556'  ='Severe', # Non-ARDS
                                'covid_557'  ='Severe', # ARDS
                                'covid_558'  ='Severe', # ARDS
                                'covid_559'  ='Severe', # ARDS
                                'covid_560'  ='Severe', # ARDS
                                'covid_561'  ='Severe', # Non-ARDS
                                'HIP002' = 'Healthy Control',
                                'HIP015' = 'Healthy Control',
                                'HIP023' = 'Healthy Control',
                                'HIP043' = 'Healthy Control',
                                'HIP044' = 'Healthy Control',
                                'HIP045' = 'Healthy Control',
                                
                                'C19-CB-0001' = 'Mild',
                                'C19-CB-0002' = 'Mild',
                                'C19-CB-0003' = 'Mild',
                                'C19-CB-0204' = 'Mild',
                                'C19-CB-0214' = 'Mild',
                                'C19-CB-0008' = 'Severe',
                                'C19-CB-0009' = 'Severe',
                                'C19-CB-0011' = 'Severe',
                                'C19-CB-0012' = 'Severe',
                                'C19-CB-0013' = 'Severe',
                                'C19-CB-0016' = 'Severe',
                                'C19-CB-0020' = 'Severe',
                                'C19-CB-0021' = 'Severe',
                                'C19-CB-0198' = 'Severe',
                                'C19-CB-0199' = 'Severe',
                                
                                'BAL009' =	'Mild',
                                'BAL012' =	'Severe',
                                'BAL013' =	'Severe',
                                'BAL014' =	'Severe',
                                'BAL015' =	'Severe',
                                'BAL016' =	'Severe',
                                'BAL022' =	'Severe',
                                'BAL023' =	'Severe',
                                'BAL024' =	'Severe',
                                'BAL025' =	'Severe',
                                'BAL032' =	'Severe',
                                'BAL035' =	'Severe',
                                'BAL039' =	'Severe',
                                'BAL040' =	'Severe',
                                'BAL027' =	'Severe',
                                'BAL031' =	'Severe',
                                'BAL037' =	'Mild',
                                'BAL021' =	'Severe',
                                'BAL026' =	'Severe',
                                'BAL033' =	'Severe',
                                'BAL034' =	'Severe',
                                'BAL020' =	'Severe',
                                
                                'GSM4593897_sample_10_UMI_counts.csv'= 'Healthy Control',
                                'GSM4593896_sample_9_UMI_counts.csv' = 'Healthy Control',
                                'GSM4593895_sample_8_UMI_counts.csv' = 'Healthy Control',
                                'GSM4593894_sample_7_UMI_counts.csv' = 'Healthy Control',
                                'GSM4593893_sample_6_UMI_counts.csv' = 'Healthy Control',
                                'GSM4593892_sample_5_UMI_counts.csv' = 'Healthy Control',
                                'GSM4593891_sample_4_UMI_counts.csv' = 'Healthy Control',
                                'GSM4593890_sample_3_UMI_counts.csv' = 'Healthy Control',
                                'GSM4593889_sample_2_UMI_counts.csv' = 'Healthy Control',
                                'GSM4593888_sample_1_UMI_counts.csv' = 'Healthy Control')

immune.combined[["group"]] <- immune.combined@active.ident

# Add compartment information
Idents(object=immune.combined) <- "study"

immune.combined <- RenameIdents(object = immune.combined,
                                'Liao'='BALF',
                                'Wauters'='BALF',
                                'Tabib'='BALF',
                                'Mould'='BALF',
                                'Lee'='PBMC',
                                'Arun'='PBMC',
                                'Schulte'='PBMC',
                                'Wilk'='PBMC')

immune.combined[["compartment"]] <- immune.combined@active.ident

# Add cohort information
Idents(object=immune.combined) <- "study"

immune.combined <- RenameIdents(object = immune.combined,
                                'Liao'='Liao',
                                'Wauters'='Wauters',
                                'Tabib'='BALFxControl',
                                'Mould'='BALFxControl',
                                'Lee'='Lee',
                                'Arun'='Arun',
                                'Schulte'='Schulte',
                                'Wilk'='Wilk')

immune.combined[["cohort"]] <- immune.combined@active.ident

# Visualize UMAP according to severity group
Idents(object=immune.combined) <- "cohort"
pdf(file="harmony_umap_group_split.pdf",width=12,height=6)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, split.by = "group",label=FALSE,raster=FALSE)
dev.off()

# Visualize UMAP according to cohort
pdf(file="harmony_umap_cohort.pdf",width=7,height=6)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "cohort",label=FALSE,raster=FALSE)
dev.off()

# Visualize UMAP according to BALF vs. PBMC
pdf(file="harmony_umap_compartment.pdf",width=6,height=6)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "compartment",label=FALSE,raster=FALSE)
dev.off()

pdf(file="harmony_umap_compartment_split.pdf",width=12,height=6)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, split.by = "compartment",label=FALSE,raster=FALSE)
dev.off()

########################### ADDITIONAL QUALITY CONTROL TO PREP FOR DE ###########################

immune.combined[['percent.hemo']] <- PercentageFeatureSet(immune.combined, pattern = "^HB[^(P)]", assay = 'RNA')
immune.combined[['percent.ppbp']] <- PercentageFeatureSet(immune.combined, pattern = "^PPBP$", assay = 'RNA')

Idents(object = immune.combined) <- "seurat_clusters"

# Plot original
p0 <- VlnPlot(immune.combined,pt.size =0.5,features="percent.hemo")
pdf(file = "pct.hemo_og.pdf")
plot(p0)
dev.off()

p0 <- VlnPlot(immune.combined,pt.size =0.5,features="percent.ppbp")
pdf(file = "pct.ppbp_og.pdf")
plot(p0)
dev.off()

# Filter cells by hemo
immune.combined <- subset(x = immune.combined, subset = percent.hemo < 5)
# Filter cells by PPBP
immune.combined <- subset(x = immune.combined, subset = percent.ppbp < 1)

# Plot again
p0 <- VlnPlot(immune.combined,pt.size =0.5,features="percent.hemo")
pdf(file = "pct.hemo_rem.pdf")
plot(p0)
dev.off()

p0 <- VlnPlot(immune.combined,pt.size =0.5,features="percent.ppbp")
pdf(file = "pct.ppbp_rem.pdf")
plot(p0)
dev.off()

########################### SAVE FINAL OBJECT ###########################

saveRDS(immune.combined,"main_fine_labeled.rds")

print("Script finished!")