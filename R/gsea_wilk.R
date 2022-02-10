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
DE_dir <- paste0(wd,"/Differential_Expression/DE_Results[1]")

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



#Function for waterfall plots

waterfall_plot <- function (fgsea_results, graph_title) {
  fgsea_results %>%
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 12,face="bold"),
          plot.title = element_text(hjust = 1))
}


# Function to perform GSEA
run_GSEA <- function(dataset,de_celltypes,groups1,groups2){
  
  # select patient group for ident1 
  for (group1 in  groups1){
    # select patient group for ident2
    for (group2 in groups2){
      # if they are not the same make a comparison
      if(group1 != group2){
        for (item in de_celltypes) {
          setwd(DE_dir)
          print(item)
          raw_list <- read.csv(file = paste0(dataset,group1,"_",group2,"_",item,"_DE_nCount_RNA.csv"),header = T, row.names = 1)
          
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
          fwrite(fgsea_results,file = paste0("GSEA Results/",dataset,item," ",group1," ",group2," fgsea.tsv"),sep="\t",sep2=c(""," ",""))
          
          p0 <- waterfall_plot(fgsea_results, paste0("Pathways Enriched in ", dataset,item ," ",group1," ",group2))
          p0
          
          pdf(file = paste0("GSEA Waterfalls/",dataset,item,"_",group1,"_",group2,"_WATERFALL.pdf"),height=9,width=7)
          plot(p0)
          dev.off()
        }
      }
    }
  }}



########################## WILK ##########################

dataset_wilk <- "Wilk_" # Enter the dataset header
de_celltypes_wilk <- c("CD14 Mono",
                       "CD16 Mono",
                       "CD4 Naive T",
                       "CD4 Treg",
                       "CD8 Effector T",
                       "CD8 Naive T",
                       # "cDC1",
                       # "cDC2",
                       "Inflammatory MP",
                       "Intermediate Mono",
                       "Mature B",
                       "Neutrophil",
                       "NK",
                       "pDC",
                       "Plasma cell",
                       "Proliferating T")

groups1_wilk <- c("Severe") # Enter the test severity categories here (FC=group1/group2)
groups2_wilk <- c("Healthy Control") # Enter the reference severity categories here

run_GSEA(dataset_wilk,de_celltypes_wilk,groups1_wilk,groups2_wilk)

warnings()
