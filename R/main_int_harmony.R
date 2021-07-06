library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)
library(Matrix)
library(DoubletFinder)
set.seed(10097)

print("Harmony all datasets, doublets removed, 0.05% + 5 PCs") # This uses 22 PCs!
print("Check all outputs are the same as before")

wd <- getwd()

########################## Arun PBMC ##########################

data_dir<-paste0(wd,"/arun")
setwd(data_dir)

arun.list <- list()
samplenames <- c("cov1","cov2","cov3","cov4","cov7","cov8","cov9","cov10","cov11","cov12","cov17","cov18")

for (sample in samplenames){
  
  # Loading in CITE-seq data (multimodal)
  data <- Read10X(sample)
  sample.tmp.seurat = CreateSeuratObject(counts = data$`Gene Expression`,min.cells = 3,min.genes = 200) # Store the gene expression data
  sample.tmp.seurat[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`) # Store the Ab data
  
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample.tmp.seurat[['percent.18S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA18S5", assay = 'RNA')
  sample.tmp.seurat[['percent.28S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA28S5", assay = 'RNA')
  
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.18S <20 & percent.28S <20 & percent.mito < 10)
  
  sample.tmp.seurat@meta.data[["orig.ident"]] <- sample
  
  arun.list <-append(arun.list, sample.tmp.seurat)
}
names(arun.list) <- samplenames
rm(sample.tmp.seurat,data)

print("Arun loaded!")


########################## LEE PBMC ##########################

# Set directory
data_dir <- paste0(wd,"/lee")
setwd(data_dir)

# Read in Lee data
rawdata <- Read10X("filtered_gene_bc_matrices")
# Numerical suffix of each barcode represents the patient information
data <- CreateSeuratObject(rawdata,min.cells = 3,min.genes = 200, names.delim = "-",names.field = 2)
rm(rawdata)

# Add donor names to metadata
Idents(object=data) <- "orig.ident"
# Note: the numbers for each sample are found directly in the GEO
data <- RenameIdents(object = data,
                     '1'  = 'nCoV 1 scRNA-seq',
                     '2'  = 'nCoV 2 scRNA-seq',
                     '3'  = 'Flu 1 scRNA-seq',
                     '4'  = 'Flu 2 scRNA-seq',
                     '5'  = 'Normal 1 scRNA-seq',
                     '6'  = 'Flu 3 scRNA-seq',
                     '7'  = 'Flu 4 scRNA-seq',
                     '8'  = 'Flu 5 scRNA-seq',
                     '9'  = 'nCoV 3 scRNA-seq',
                     '10' = 'nCoV 4 scRNA-seq',
                     '11' = 'nCoV 5 scRNA-seq',
                     '12' = 'nCoV 6 scRNA-seq',
                     '13' = 'Normal 2 scRNA-seq',
                     '14' = 'Normal 3 scRNA-seq',
                     '15' = 'nCoV 7 scRNA-seq',
                     '16' = 'nCoV 8 scRNA-seq',
                     '17' = 'nCoV 9 scRNA-seq',
                     '18' = 'nCoV 10 scRNA-seq',
                     '19' = 'Normal 4 scRNA-seq',
                     '20' = 'nCoV 11 scRNA-seq')

data[["donor"]] <- data@active.ident
data[["orig.ident"]] <- data@active.ident # Rename the original idents

obj.list <- SplitObject(data, split.by = "donor")

# Preprocess
lee.list = list()
for(sample.tmp.seurat in obj.list){
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample.tmp.seurat[['percent.18S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA18S5", assay = 'RNA')
  sample.tmp.seurat[['percent.28S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA28S5", assay = 'RNA')
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.18S <20 & percent.28S <20 & percent.mito < 10)
  lee.list <-append(lee.list, sample.tmp.seurat)
}
names(lee.list) <- names(obj.list)
rm(obj.list,data,sample.tmp.seurat)

# Remove the flu samples from the dataset
lee.list[['Flu 1 scRNA-seq']] <- NULL
lee.list[['Flu 2 scRNA-seq']] <- NULL
lee.list[['Flu 3 scRNA-seq']] <- NULL
lee.list[['Flu 4 scRNA-seq']] <- NULL
lee.list[['Flu 5 scRNA-seq']] <- NULL
# Remove the asymptomatic patient
lee.list[['nCoV 11 scRNA-seq']] <- NULL
# Remove repeated later timepoints
lee.list[['nCoV 4 scRNA-seq']] <- NULL # Remove later timepoint for patient C3
lee.list[['nCoV 8 scRNA-seq']] <- NULL # Remove later timepoint for patient C6
lee.list[['nCoV 10 scRNA-seq']] <- NULL # Remove later timepoint for patient C7

print("Lee loaded!")


# ########################## WILK PBMC ##########################

# Set directory
data_dir<-paste0(wd,"/wilk")
setwd(data_dir)

#Read in the Blish Seurat object
raw <- readRDS("blish_covid.seu.rds")
DefaultAssay(raw) <- "RNA"
UpdateSeuratObject(raw)

#Split into a list of subsetted objects
obj.list <- SplitObject(raw, split.by = "orig.ident")

# Preprocess
wilk.list = list()
for(sample.tmp.seurat in obj.list){
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample.tmp.seurat[['percent.18S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA18S5", assay = 'RNA')
  sample.tmp.seurat[['percent.28S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA28S5", assay = 'RNA')
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.18S <20 & percent.28S <20 & percent.mito < 10)
  wilk.list <-append(wilk.list, sample.tmp.seurat)
}

names(wilk.list) <- names(obj.list)

wilk.list["covid_555_2"] <- NULL # Remove the later timepoint for patient C1 (C1B)

rm(obj.list,raw,sample.tmp.seurat)

print("Wilk loaded!")

########################## SCHULTE PBMC ##########################

# Set directory
data_dir<-paste0(wd,"/schulte")
setwd(data_dir)

# Read in the Schulte Seurat object
raw <- readRDS("seurat_COVID19_PBMC_cohort1_10x.rds")
DefaultAssay(raw) <- "RNA"
UpdateSeuratObject(raw)

raw@meta.data[["orig.ident"]] <- raw@meta.data[["donor"]] # orig.ident wasn't named uniquely, stash donor names as orig.ident

#Split into a list of subsetted objects
obj.list <- SplitObject(raw, split.by = "orig.ident")

# Delete frozen PBMCs only COVID donors
obj.list["C19-CB-0005"] <- NULL        
obj.list["C19-CB-0052"] <- NULL        
obj.list["C19-CB-0053"] <- NULL 

# Get rid of the controls that came from Miguel's paper
obj.list["P18F"] <- NULL        
obj.list["P17H"] <- NULL        
obj.list["P20H"] <- NULL        
obj.list["P15F"] <- NULL        
obj.list["P08H"] <- NULL        
obj.list["P13H"] <- NULL        
obj.list["P07H"] <- NULL        
obj.list["P06F"] <- NULL       
obj.list["P04H"] <- NULL
obj.list["C2P01H"] <- NULL      
obj.list["P09H"] <- NULL        
obj.list["P02H"] <- NULL        
obj.list["C2P05F"] <- NULL      
obj.list["C2P07H"] <- NULL      
obj.list["C2P13F"] <- NULL      
obj.list["C2P16H"] <- NULL      
obj.list["C2P10H"] <- NULL      
obj.list["C2P19H"] <- NULL      
obj.list["C2P15H"] <- NULL

# Get rid of the standard 10x runs
obj.list['one_k_v3'] <- NULL
obj.list['Five_k_v3'] <- NULL
obj.list['Ten_k_v3'] <- NULL

# Preprocess
schulte.list = list()
for(sample.tmp.seurat in obj.list){
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample.tmp.seurat[['percent.18S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA18S5", assay = 'RNA')
  sample.tmp.seurat[['percent.28S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA28S5", assay = 'RNA')
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.18S <20 & percent.28S <20 & percent.mito < 10)
  # The original object already was filtered for genes expressed in > 5 cells
  schulte.list <-append(schulte.list, sample.tmp.seurat)
}

names(schulte.list) <- names(obj.list)

# Take only the sample from the earliest timepoint
temp <- schulte.list[["C19-CB-0008"]]
Idents(object=temp) <- "days_after_onset"
schulte.list[["C19-CB-0008"]] <- subset(temp,idents='13')

temp <- schulte.list[["C19-CB-0009"]]
Idents(object=temp) <- "days_after_onset"
schulte.list[["C19-CB-0009"]] <- subset(temp,idents='9')

temp <- schulte.list[["C19-CB-0012"]]
Idents(object=temp) <- "days_after_onset"
schulte.list[["C19-CB-0012"]] <- subset(temp,idents='9')

temp <- schulte.list[["C19-CB-0013"]]
Idents(object=temp) <- "days_after_onset"
schulte.list[["C19-CB-0013"]] <- subset(temp,idents='8')

rm(obj.list,raw,sample.tmp.seurat,temp)


print("Schulte loaded!")


# ########################## LIAO BALF ##########################

# Set directory
data_dir<-paste0(wd,"/liao")
setwd(data_dir)

#Read in the Zhang BALF .h5 files downloaded from the GEO
samples_h5 = read.delim2("meta_h5.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
samples_tsv = read.delim2("meta_tsv.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")

# sample_i stores the sample info
liao.list = list()
for(sample_s in samples_h5$sample){
  print(sample_s)
  sample_i = samples_h5 %>% dplyr::filter(.,sample == sample_s)
  filename = paste(sample_s,"_filtered_feature_bc_matrix.h5",sep="")
  print(filename)
  sample.tmp = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  sample.tmp.seurat <- CreateSeuratObject(counts = sample.tmp, min.cells = 3, min.features = 200,project = sample_s)
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample.tmp.seurat[['percent.18S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA18S5", assay = 'RNA')
  sample.tmp.seurat[['percent.28S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA28S5", assay = 'RNA')
  sample_i$nFeature_RNA_low = as.numeric(sample_i$nFeature_RNA_low)
  sample_i$nFeature_RNA_high = as.numeric(sample_i$nFeature_RNA_high)
  sample_i$nCount_RNA = as.numeric(sample_i$nCount_RNA)
  sample_i$percent.mito = as.numeric(sample_i$percent.mito)
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.18S <20 & percent.28S <20 & percent.mito < 10)
  liao.list[sample_s] = sample.tmp.seurat
}

for(sample_s in samples_tsv$sample){
  print(sample_s)
  sample_i = samples_tsv %>% dplyr::filter(.,sample == sample_s)
  sample.tmp = Read10X(sample_s)
  sample.tmp.seurat <- CreateSeuratObject(counts = sample.tmp, min.cells = 3, min.features = 200,project = sample_s)
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample.tmp.seurat[['percent.18S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA18S5", assay = 'RNA')
  sample.tmp.seurat[['percent.28S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA28S5", assay = 'RNA')
  sample_i$nFeature_RNA_low = as.numeric(sample_i$nFeature_RNA_low)
  sample_i$nFeature_RNA_high = as.numeric(sample_i$nFeature_RNA_high)
  sample_i$nCount_RNA = as.numeric(sample_i$nCount_RNA)
  sample_i$percent.mito = as.numeric(sample_i$percent.mito)
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.18S <20 & percent.28S <20 & percent.mito < 10)
  
  liao.list[sample_s] = sample.tmp.seurat
}
print("Liao loaded!")

########################## WAUTERS BALF ##########################

# Set directory
data_dir<-paste0(wd,"/wauters")
setwd(data_dir)

wauters_dgc <- readRDS("wauters_Allcells.counts.rds")
wauters_orig <- CreateSeuratObject(counts = wauters_dgc, min.cells = 3, min.features = 200, project = "wauters")
rm(wauters_dgc)

wauters_orig[["donor"]] <- wauters_orig@active.ident

metadata <- read.csv("wauters_Allcells.meta.data.csv")
metadata<-unique(metadata[,c(2:4)])
covid_dnrs <- metadata[metadata$Disease=="COVID19",]$PatientNumber

wauters = subset(wauters_orig,idents = covid_dnrs) # Subset for just the COVID donors
rm(wauters_orig)
rm(metadata)

DefaultAssay(wauters) <- "RNA"
UpdateSeuratObject(wauters)

obj.list <- SplitObject(wauters, split.by = "orig.ident")

print("Wauters object names")
print(names(obj.list))

# Preprocess
wauters.list = list()
for(sample.tmp.seurat in obj.list){
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample.tmp.seurat[['percent.18S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA18S5", assay = 'RNA')
  sample.tmp.seurat[['percent.28S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA28S5", assay = 'RNA')
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.18S <20 & percent.28S <20 & percent.mito < 10)
  wauters.list <-append(wauters.list, sample.tmp.seurat)
}
names(wauters.list) <- names(obj.list)

rm(obj.list,sample.tmp.seurat)
print("Wauters loaded!")

########################## BALF CROSS-CONTROLS ##########################
controls.list = list()

filenames <- c("GSM4593888_sample_1_UMI_counts.csv",
               "GSM4593889_sample_2_UMI_counts.csv",
               "GSM4593890_sample_3_UMI_counts.csv",
               "GSM4593891_sample_4_UMI_counts.csv",
               "GSM4593892_sample_5_UMI_counts.csv",
               "GSM4593893_sample_6_UMI_counts.csv",
               "GSM4593894_sample_7_UMI_counts.csv",
               "GSM4593895_sample_8_UMI_counts.csv",
               "GSM4593896_sample_9_UMI_counts.csv",
               "GSM4593897_sample_10_UMI_counts.csv")

for(filename in filenames){
  print(filename)
  sample.tmp = read.csv(filename,row.names=1)
  sample.tmp.seurat <- CreateSeuratObject(counts = sample.tmp, min.cells = 3, min.features = 200,project = filename)
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample.tmp.seurat[['percent.18S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA18S5", assay = 'RNA')
  sample.tmp.seurat[['percent.28S']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^RNA28S5", assay = 'RNA')
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1000 & percent.18S <20 & percent.28S <20 & percent.mito < 10)
  controls.list[filename] = sample.tmp.seurat
}
names(controls.list) <- filenames

rm(sample.tmp,sample.tmp.seurat,wauters)

print("Controls loaded!")


########################## CONCATENATE THE LISTS ###########################
# Concatenate the lists!
full.list <- c(liao.list,lee.list,arun.list,wilk.list,wauters.list,controls.list,schulte.list)

# Check that each object has a name and check number of cells
print(sapply(full.list, ncol))

# Clean up the workspace
rm(liao.list,lee.list,arun.list,wilk.list,wauters.list,controls.list,schulte.list)


########################## DOUBLET REMOVAL ###########################

for(i in 1:length(full.list)){
  
  # Get sample name (filename)
  minorList <- full.list[i]
  filename <- names(minorList)
  
  print(filename)
  
  # IMPORTANT: GET TWO SEURAT OBJECTS
  original.seurat <- minorList[[1]]
  doublet.seurat <- minorList[[1]]
  
  ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
  doublet.seurat <- NormalizeData(doublet.seurat)
  doublet.seurat <- ScaleData(doublet.seurat)
  doublet.seurat <- FindVariableFeatures(doublet.seurat, selection.method = "vst", nfeatures = 2000)
  doublet.seurat <- RunPCA(doublet.seurat)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.doublet <- paramSweep_v3(doublet.seurat, PCs = 1:10, sct = FALSE)
  sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = FALSE)
  bcmvn.doublet <- find.pK(sweep.stats.doublet)
  pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$BCmetric)]; pK <- as.numeric(levels(pK))[pK]
  
  EDFR=0.075 # Estimated doublet formation rate
  nExp_poi <- round(EDFR*nrow(doublet.seurat@meta.data)) 
  
  doublet.seurat <- doubletFinder_v3(doublet.seurat, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_poi)
  
  attribute <- paste('pANN', 0.25, pK, nExp_poi, sep = '_')
  class <- paste('DF.classifications', 0.25, pK, nExp_poi, sep = '_')
  doublet.seurat@meta.data[['DoubletClass']] <- doublet.seurat@meta.data[[class]]
  
  score <- doublet.seurat@meta.data[[attribute]]
  
  # calculate threshold based on identification rate equal to the true doublet rate
  threshold <- sort(score, decreasing = TRUE)[nExp_poi]
  n_removed <- ncol(subset(doublet.seurat, DoubletClass=="Doublet"))
  print(threshold)
  print(n_removed)
  
  # Put doublet annotations into "original.seurat"
  original.seurat@meta.data[['DoubletClass']] <- doublet.seurat@meta.data[['DoubletClass']]
  original.seurat@meta.data[['DoubletScore']] <- doublet.seurat@meta.data[[attribute]]
  
  # REMOVE DOUBLETS in original.seurat, Save original.seurat back into full.list
  full.list[[i]] <- subset(original.seurat, DoubletClass=="Singlet") # ESSENTAIL step for removing doublets
  
}


########################## PREPROCESS AND RUN PCA ###########################

# Return working directory to root
setwd(wd)

# Conduct test to choose the number of PCs to include
pca.test<- merge(full.list[[1]], y = full.list[2:length(full.list)],project = "Harmony_alldata") %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))

# Remove mitochondrial, rRNA, and ribo genes from the variable genes list
var.genes <- pca.test@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', pca.test@assays$RNA@var.features)]
pca.test <- RunPCA(pca.test, features = var.genes, npcs = 100, verbose = FALSE) # Start by calculating 100 PCs

# Choose number of PCs 
pca.test <- ProjectDim(object = pca.test)
pdf(file="test_harmony_elbowplot.pdf")
ElbowPlot(object = pca.test,ndims = 100)
dev.off()

PCTthresh <- 0.05 # Set threshold here for determining the number of PCs to use
n_beyond <- 5 # Set the number of PCs to use beyond the elbow
pct <- pca.test[["pca"]]@stdev / sum(pca.test[["pca"]]@stdev) * 100
nPCs <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > PCTthresh), decreasing = T)[1] + 1 + n_beyond

cat("/nNumber of PCs used for Harmony: ")
print(nPCs)

# Preprocess the immune.combined object to run the actual PCA
immune.combined <- merge(full.list[[1]], y = full.list[2:length(full.list)],project = "Harmony_alldata") %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE,vars.to.regress = c("percent.mito","nCount_RNA"))

# Print the mitochondrial, rRNA, and ribo genes removed
removed.genes <- immune.combined@assays$RNA@var.features[grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', immune.combined@assays$RNA@var.features)]
print(paste("Removed genes: ",removed.genes))

# Remove the genes from the immune.combined object, run the actual PCA using the number of PCs calculated in the test
var.genes <- immune.combined@assays$RNA@var.features[!grepl('^MT-|^RPS|^RPL|^RNA18S|^RNA28S', immune.combined@assays$RNA@var.features)]
immune.combined <- RunPCA(immune.combined, features = var.genes, npcs = nPCs, verbose = FALSE)

print("Scale and PCA successful!")


########################## ADD STUDY METADATA ###########################

Idents(object=immune.combined) <- "orig.ident"

immune.combined <- RenameIdents(object = immune.combined,
                     'C51'='Liao',                                
                     'C52'='Liao',
                     'C100'='Liao',                              
                     'C141'='Liao',
                     'C142'='Liao',                              
                     'C144'='Liao',
                     'C145'='Liao',                              
                     'C143'='Liao',
                     'C146'='Liao',                              
                     'C148'='Liao',
                     'C149'='Liao',                              
                     'C152'='Liao',
                     'GSM3660650'='Tabib',                  
                     'nCoV 1 scRNA-seq'='Lee',
                     'nCoV 2 scRNA-seq'='Lee',               
                     'Normal 1 scRNA-seq'='Lee',
                     'nCoV 3 scRNA-seq'='Lee',                  
                     'nCoV 5 scRNA-seq'='Lee',                    
                     'nCoV 6 scRNA-seq'='Lee',
                     'Normal 2 scRNA-seq'='Lee',                
                     'Normal 3 scRNA-seq'='Lee',
                     'nCoV 7 scRNA-seq'='Lee',                 
                     'nCoV 9 scRNA-seq'='Lee',               
                     'Normal 4 scRNA-seq'='Lee',                 
                     'cov1'='Arun',                              
                     'cov2'='Arun',
                     'cov3'='Arun',                              
                     'cov4'='Arun',
                     'cov7'='Arun',                              
                     'cov8'='Arun',
                     'cov9'='Arun',                             
                     'cov10'='Arun',
                     'cov11'='Arun',                              
                     'cov12'='Arun',
                     'cov17'='Arun',                              
                     'cov18'='Arun',
                     'covid_555_1'='Wilk',                        
                     'covid_556'='Wilk',                          
                     'covid_557'='Wilk',
                     'covid_558'='Wilk',                          
                     'covid_559'='Wilk',
                     'covid_560'='Wilk',                          
                     'covid_561'='Wilk',
                     'HIP002'='Wilk',                            
                     'HIP015'='Wilk',
                     'HIP023'='Wilk',                            
                     'HIP043'='Wilk',
                     'HIP044'='Wilk',                            
                     'HIP045'='Wilk',
                     'C19-CB-0001'='Schulte',
                     'C19-CB-0002'='Schulte',
                     'C19-CB-0003'='Schulte',
                     'C19-CB-0204'='Schulte',
                     'C19-CB-0214'='Schulte',
                     'C19-CB-0008'='Schulte',
                     'C19-CB-0009'='Schulte',
                     'C19-CB-0011'='Schulte',
                     'C19-CB-0012'='Schulte',
                     'C19-CB-0013'='Schulte',
                     'C19-CB-0016'='Schulte',
                     'C19-CB-0020'='Schulte',
                     'C19-CB-0021'='Schulte',
                     'C19-CB-0198'='Schulte',
                     'C19-CB-0199'='Schulte',
                     'BAL009'='Wauters',                            
                     'BAL012'='Wauters',
                     'BAL013'='Wauters',                            
                     'BAL014'='Wauters',
                     'BAL015'='Wauters',                            
                     'BAL016'='Wauters',
                     'BAL022'='Wauters',                            
                     'BAL023'='Wauters',
                     'BAL024'='Wauters',                            
                     'BAL025'='Wauters',
                     'BAL032'='Wauters',                            
                     'BAL035'='Wauters',
                     'BAL039'='Wauters',                            
                     'BAL040'='Wauters',
                     'BAL027'='Wauters',                            
                     'BAL031'='Wauters',
                     'BAL037'='Wauters',                            
                     'BAL021'='Wauters',
                     'BAL026'='Wauters',                            
                     'BAL033'='Wauters',
                     'BAL034'='Wauters',                            
                     'BAL020'='Wauters',
                     'GSM4593897_sample_10_UMI_counts.csv'='Mould', 
                     'GSM4593896_sample_9_UMI_counts.csv'='Mould',
                     'GSM4593895_sample_8_UMI_counts.csv'='Mould',
                     'GSM4593894_sample_7_UMI_counts.csv'='Mould',
                     'GSM4593893_sample_6_UMI_counts.csv'='Mould',
                     'GSM4593892_sample_5_UMI_counts.csv'='Mould',
                     'GSM4593891_sample_4_UMI_counts.csv'='Mould',
                     'GSM4593890_sample_3_UMI_counts.csv'='Mould',
                     'GSM4593889_sample_2_UMI_counts.csv'='Mould',
                     'GSM4593888_sample_1_UMI_counts.csv'='Mould')
           
immune.combined[["study"]] <- immune.combined@active.ident


########################## RUN THE HARMONY PIPELINE ###########################

# Print initial info prior to running Harmony
pdf(file="test_initial_pca.pdf")
DimPlot(object = immune.combined, reduction = "pca", pt.size = .1, group.by = "orig.ident") + NoLegend()
dev.off()

pdf(file="test_initial_violin.pdf")
VlnPlot(object = immune.combined, features = "PC_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
dev.off()

# RUN HARMONY, make sure to specify what variable(s) to regress out here
# immune.combined <- immune.combined %>% RunHarmony(group.by.vars=c("orig.ident","study"),theta=c(2,0), plot_convergence = FALSE)
immune.combined <- immune.combined %>% RunHarmony("orig.ident", plot_convergence = FALSE)

# Let's make sure that the datasets are well integrated in the first 2 dimensions after Harmony.
pdf(file="test_harmony_pca.pdf")
DimPlot(object = immune.combined, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + NoLegend()
dev.off()

pdf(file="test_harmony_violins.pdf")
VlnPlot(object = immune.combined, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
dev.off()

# Clustering
immune.combined <- immune.combined %>% 
  RunUMAP(reduction = "harmony", dims = 1:nPCs) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:nPCs) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

print("\nNumber of donors: ")
print(length(unique(immune.combined@meta.data[["orig.ident"]])))

print("\nTotal number of cells: ")
print(ncol(immune.combined@assays[["RNA"]]))

# Visualization
pdf(file="test_harmony_umap_clusters.pdf")
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = .1)
dev.off()

pdf(file="test_harmony_umap_study.pdf",width=6,height=6)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "study",label=FALSE)
dev.off()

pdf(file="test_harmony_umap_donors_showall.pdf",width=15,height=6)
DimPlot(object = immune.combined, reduction = "umap", pt.size = .1, group.by = "orig.ident",label=FALSE)
dev.off()

saveRDS(immune.combined,"test_main_harmony_integrated.rds")

########################## MARKER DOTPLOTS AND HEATMAPS ###########################

Idents(object=immune.combined) <- "seurat_clusters"
DefaultAssay(immune.combined) <- "RNA" 

# Dotplots for cell type annotation
features=rev(c("EMP2","TPPP3","KRT18","KRT15","SCGB1A1","TPSB2","TPSAB1","TYMS","MKI67","IGHA1","IGHG4","IGHM",
               "IGHD","CD22","CD19","MS4A1","CD79A","CXCR6","GZMK","KLRB1","EOMES","KLRF1","GZMB","GNLY","KLRD1",
               "CD69","FAS","FOXP3","CTLA4","IL2RA","CD44","IL7R","IL2RG","TCF7","CCR7","SELL","CD8A","CD4","CD3D",
               "CXCR2","FCGR3B","HCAR3","PI3","ELANE","LILRA4","IL3RA","THBD","CD1C","CLEC9A","XCR1","FABP4","SPP1",
               "FCGR3A","CD14","S100A8","FCN1","CD68","PPBP")) 

p4 <- DotPlot(immune.combined, features = features,dot.scale=8) # + RotatedAxis()
pdf(file = "test_main_dotplot.pdf",height=10,width=30)
plot(p4)
p4
dev.off()

print("Script finished!")