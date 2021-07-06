# For MAST categorical variables, see https://github.com/satijalab/seurat/issues/3884

require("knitr")
library(stringr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(future)
library(tidyverse)
library(data.table)
library(EnhancedVolcano)
set.seed(10097)

wd <- getwd()
res_dir<-paste0(wd,"/Differential_Expression")
print(res_dir)
opts_knit$set(root.dir = res_dir)

# Create a new directory for DE Results (won't overwrite an existing directory)
dir.create(file.path(res_dir, "DE_Results[5]"), showWarnings = TRUE)
dir.create(file.path(res_dir, "DE_Volcanos[5]"), showWarnings = TRUE)


#Define function to perform DE

get_DE_genes <- function (seuratobj, groups, covars) {
  seuratobj$donor <- droplevels(x = seuratobj$donor) # Drop the excess donor levels if using "donor" as a MAST covariate
  
  DefaultAssay(seuratobj) <- "RNA"
  
  # select patient severity group
  for (group in  groups){
    
    # Select the appropriate studies for the comparison
    if (group=="Healthy Control"){
      balf_studies <- c("Liao","Tabib","Mould")
      pbmc_studies <- c("Wilk","Arun","Lee")
      
      de_celltypes <- c(
        "CD14 Mono",
        "CD16 Mono",
        "CD4 Naive T",
        "CD4 Treg",
        "CD8 Effector T",
        "CD8 Naive T",
        "cDC1",
        "cDC2",
        "Inflammatory MP",
        # "Intermediate Mono",
        "Mature B",
        "Neutrophil",
        "NK",
        "pDC",
        "Plasma cell",
        "Proliferating T"
      )
    }
    else if (group=="Mild"){
      balf_studies <- c("Liao","Wauters")
      pbmc_studies <- c("Arun","Lee","Schulte")
      
      de_celltypes <- c(
        "CD14 Mono",
        "CD16 Mono",
        "CD4 Naive T",
        "CD4 Treg",
        "CD8 Effector T",
        "CD8 Naive T",
        "cDC1",
        "cDC2",
        "Inflammatory MP",
        # # "Intermediate Mono",
        "Mature B",
        "Neutrophil",
        "NK",
        # # "pDC",
        # "Plasma cell",
        "Proliferating T"
      )
    }
    else if (group=="Severe"){
      balf_studies <- c("Liao","Wauters")
      pbmc_studies <- c("Wilk","Arun","Lee","Schulte")
      
      de_celltypes <- c(
        "CD14 Mono",
        "CD16 Mono",
        "CD4 Naive T",
        "CD4 Treg",
        "CD8 Effector T",
        "CD8 Naive T",
        "cDC1",
        "cDC2",
        "Inflammatory MP",
        # "Intermediate Mono",
        "Mature B",
        "Neutrophil",
        "NK",
        "pDC",
        "Plasma cell",
        "Proliferating T"
      )
    }
    
    
    # loop through each of the cell types we want DE for
    for (item in de_celltypes) {
      print(item)

      # Subset out the celltype of interest
      Idents(object=seuratobj) <- "celltype_fine"
      cell_subset <- subset(seuratobj, idents = item)
      
      # Subset out the severity group of interest
      Idents(object=cell_subset) <- "group"
      de_subset <- subset(cell_subset, idents = group)
      
      # Set identities to compartment (BALF/PBMC)
      de_subset <- SetIdent(de_subset, value = "compartment") 

      de_results <- FindMarkers(de_subset, "BALF", "PBMC", logfc.threshold = 0, min.pct=0, min.cells.group=10, test.use="MAST",latent.vars=covars)

      write.csv(de_results[order(de_results$avg_log2FC, decreasing = T),],
                file = paste0("Differential_Expression/DE_Results[5]/",group,"_",item,"_DE_",paste(covars, collapse='_'),".csv"))

      volcano_plot(de_results,group,item,covars)
        
      
      }
    }
  }


volcano_plot <- function(markers,group,item,covars){
  
  p1 <- EnhancedVolcano(markers,lab = rownames(markers),title =paste0(group,"_",item)
                        ,x= 'avg_log2FC',
                        y = 'p_val_adj', 
                        xlab = bquote(~Log[2]~ '(fold change)'),
                        pCutoff = 10e-3,
                        FCcutoff = 1,
                        gridlines.major = FALSE, 
                        gridlines.minor = FALSE,
                        #pointSize = 6.0,
                        subtitle = NULL, 
                        #axisLabSize = 40,
                        #drawConnectors = TRUE,
                        #labSize = 7, 
                        #legendLabSize = 30,
                        col = c("grey","grey","grey","red3"),colAlpha=1)
  p1
  pdf(file = paste0("Differential_Expression/DE_Volcanos[5]/",group,"_",item,"_VOLCANO",paste(covars, collapse='_'),".pdf"))
  plot(p1)
  dev.off()
  
}

#------------------------------------------------------------------------------DIFFERENTIAL EXPRESSION------------------------------------------------------------------------------------------- 

fullObj <- readRDS("main_fine_labeled.rds")

# Set up multicore parallelization
plan("multicore") 
plan()
options(future.globals.maxSize = 1000 * 1024^9,seed=TRUE)


########################## DE WITHIN SAME CELL TYPE AND SEVERITY LEVEL, BETWEEN COMPARTMENTS [5] ##########################

covars=c("nCount_RNA","study")

# groups <- c("Severe","Mild","Helathy Control")
groups <- c("Mild")

print(groups)
print(covars)

# Function call
get_DE_genes(fullObj,groups,covars)

print("Script finished!")
warnings()