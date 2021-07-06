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
dir.create(file.path(res_dir, "DE_Results[3]"), showWarnings = TRUE)
dir.create(file.path(res_dir, "DE_Volcanos[3]"), showWarnings = TRUE)


#Define function to perform DE

get_DE_genes <- function (seuratobj, de_celltypes,condition1,condition2,data_type,wilk_ards,covars) {
  seuratobj$donor <- droplevels(x = seuratobj$donor) # Drop the excess donor levels if using "donor" as a MAST covariate
  groups1 <- condition1
  groups2 <- condition2
  
  Idents(object=seuratobj) <- "celltype_fine"

  DefaultAssay(seuratobj) <- "RNA"
  
  # select patient group for ident1 
  for (group1 in  groups1){
    # select patient group for ident2
    for (group2 in groups2){
      # if they are not the same make a comparison
      if(group1 != group2){
        # loop through each of the cell types you want DE for
        for (item in de_celltypes) {
          print(item)
          
          de_subset <- subset(seuratobj, idents = item)
          
          if (wilk_ards==0){
            de_subset <- SetIdent(de_subset, value = "group") # use "group" for Healthy/Mild/Severe comparisons
          }
          else if (wilk_ards==1){
            de_subset <- SetIdent(de_subset, value = "strata") # use "strata" if you want Wilk ARDS vs Non-ARDS
          }
          
          de_results <- FindMarkers(de_subset, group1, group2, logfc.threshold = 0, min.pct=0.1, min.cells.group=10, test.use="MAST",latent.vars=covars)
          
          write.csv(de_results[order(de_results$avg_log2FC, decreasing = T),],
                    file = paste0("Differential_Expression/DE_Results[3]/",data_type,"_",group1, "_",group2,"_",item,"_DE_",paste(covars, collapse='_'),".csv"))
          
          volcano_plot(de_results,data_type,group1,group2,item)
        }
      }
    }
  }
}


volcano_plot <- function(markers,data_type,group1,group2,item){
  
  p1 <- EnhancedVolcano(markers,lab = rownames(markers),title =paste0(data_type,"_",group1, "_",group2,"_",item)
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
  pdf(file = paste0("Differential_Expression/DE_Volcanos[3]/",data_type,"_",group1, "_",group2,"_",item,"_VOLCANO_",paste(covars, collapse='_'),".pdf"))
  plot(p1)
  dev.off()
  
}

#------------------------------------------------------------------------------DIFFERENTIAL EXPRESSION------------------------------------------------------------------------------------------- 

fullObj <- readRDS("main_fine_labeled.rds")

# Set up multicore parallelization
plan("multicore") 
plan()
options(future.globals.maxSize = 1000 * 1024^9,seed=TRUE)

de_celltypes <- c(
  "CD14 Mono",
  "CD16 Mono",
  "CD4 Naive T",
  "CD4 Treg",
  "CD8 Effector T",
  "CD8 Naive T",
  # "cDC1",
  "cDC2",
  "Inflammatory MP",
  # "Intermediate Mono",
  "Mature B",
  "Neutrophil",
  "NK",
  # "pDC",
  "Plasma cell",
  "Proliferating T"
)


########################## DE WITHIN DATASET SAME COMPARTMENT ##########################
cohort <- "Lee"
covars=c("nCount_RNA","donor")

print(covars)

Idents(object=fullObj) <- "cohort"
subset <- subset(fullObj, idents = cohort)

get_DE_genes(subset,de_celltypes,"Severe","Mild",cohort,0,covars)
print("Severe vs. Mild Done!")
get_DE_genes(subset,de_celltypes,"Severe","Healthy Control",cohort,0,covars)
print("Severe vs. Healthy Control Done!")
get_DE_genes(subset,de_celltypes,"Mild","Healthy Control",cohort,0,covars)
print("Mild vs. Healthy Control Done!")

warnings()
print("Script finished!")
