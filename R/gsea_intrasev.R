require("knitr")
library(stringr)
library(Seurat)
library(fgsea)
library(tidyverse)
library(future)
library(data.table)
set.seed(10097)

wd <- getwd()
res_dir<-paste0(wd,"/GSEA")
DE_dir <- paste0(wd,"/Differential_Expression/DE_Results[5]")

opts_knit$set(root.dir = res_dir)

# Create a new directory for GSEA Waterfalls (won't overwrite an existing directory)
dir.create(file.path(res_dir, "GSEA Waterfalls"), showWarnings = FALSE)
dir.create(file.path(res_dir, "GSEA Results"), showWarnings = FALSE)

# Import required gene sets
hallmark_pathway <- gmtPathways("h.all.v7.1.symbols.gmt")


# Function to format ranked list of genes for GSEA

prepare_ranked_list <- function(raw_list) {
  raw_list$Gene.name <- str_to_upper(rownames(raw_list))
  raw_list$ngLOG10p <- -log10(raw_list$p_val_adj) * sign(raw_list$avg_log2FC)
  raw_list <- raw_list[,c("Gene.name", "avg_log2FC","ngLOG10p")]
  rownames(raw_list) <- NULL
  
  # if duplicate gene names present, average the values
  if( sum(duplicated(raw_list$Gene.name)) > 0) {
    raw_list <- aggregate(.~Gene.name, FUN = mean, data = raw_list)
  }
  
  # Order the gene list (decreasing order)
  # This will not affect the result since fgsea sorts the list itself--we just do this for visualization purposes
  ranked_list <- raw_list[order(raw_list$avg_log2FC, raw_list$ngLOG10p, decreasing = T),]
  
  # omit rows with NA values
  ranked_list <- ranked_list[,c("Gene.name", "avg_log2FC")] # Note that the column you take, not the order, dictates the result
  ranked_list <- na.omit(ranked_list)
  
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  
  # check ranks are stored as a numeric vector
  print(class(ranked_list))
  
  # the names of the vector are the genes, the values are the rank metric
  print(head(ranked_list))
  
  # check ranks are from highest to lowest
  print(tail(unname(ranked_list)))
  
  return(ranked_list)
}


# Function to perform GSEA
run_GSEA <- function(condition,de_celltypes){
  
  for (item in de_celltypes) {
    setwd(DE_dir)
    print(item)
    raw_list <- read.csv(file = paste0(condition,item,"_DE_nCount_RNA_study.csv"),header = T, row.names = 1)
    
    ranked_list <- prepare_ranked_list(raw_list)
    
    set.seed(10097)
    fgsea_results <- fgseaMultilevel(pathways = hallmark_pathway,
                                     stats = ranked_list,
                                     minSize = 15,
                                     maxSize = 500,
                                     nPermSimple = 10000)
    #nperm= 1000) # Use this argument if running fgseaSimple
    
    fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()
    
    setwd(res_dir)
    fwrite(fgsea_results,file = paste0("GSEA Results/",condition,item," fgsea.tsv"),sep="\t",sep2=c(""," ",""))
    }
}



########################## RUN GSEA ##########################
celltypes_HC <- c(
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
  "Proliferating T")

celltypes_M <- c(
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
  "Proliferating T")

celltypes_S <- c(
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
  "Proliferating T")

run_GSEA("Healthy Control_",celltypes_HC)
run_GSEA("Mild_",celltypes_M)
run_GSEA("Severe_",celltypes_S)

warnings()