install.packages("enrichR")

require("knitr")
library(stringr)
library(enrichR)
library(tidyverse)
set.seed(10097)

wd <- getwd()
res_dir<-paste0(wd,"/GO/EnrichR_Results/")
DE_dir <- paste0(wd,"/Differential_Expression/DE_Plots[5]/")
dbs <- listEnrichrDbs()

opts_knit$set(root.dir = res_dir)


websiteLive <- TRUE
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)


up_down <- c("UP","DOWN")

dbs <- c(
  "GO_Molecular_Function_2018",
  "GO_Biological_Process_2018",
  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"
)

# Function to perform enrichR

enrich_fun <- function(conditions,celltype,genes,direction,dbs){
  # print(genes)
  # print(dbs)
  if (websiteLive && !is.na(genes)){
    prefix <- paste(conditions[1],conditions[2],sep = "_")
    enriched <- enrichr(genes, databases = dbs)
    
    print(enriched)
    write.table(
      enriched[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]],
      file = paste0(res_dir,"TF_",direction,"_",prefix,"_",celltype,'.txt'),
      sep = "\t",
      dec = ".",
      row.names = TRUE,
      col.names = TRUE,
      quote = FALSE
    )
    
    write.table(
      enriched[["GO_Biological_Process_2018"]],
      file = paste0(res_dir,"GO_BIO_",direction,"_",prefix,"_",celltype,".txt"),
      sep = "\t",
      dec = ".",
      row.names = TRUE,
      col.names = TRUE,
      quote = FALSE
    )
    
    write.table(
      enriched[["GO_Molecular_Function_2018"]],
      file = paste0(res_dir,"GO_Mol_",direction,"_",prefix,"_",celltype,".txt"),
      sep = "\t",
      dec = ".",
      row.names = TRUE,
      col.names = TRUE,
      quote = FALSE
    )
  }
}

# Function that takes in DE results and calls enrich_fun

enrich_preprocess <- function(conditions,full){
  condition1 <- conditions[1]
  condition2 <- conditions[2]
  
  # Truncate dataframes to recover just celltypes that match across the conditions we care about
  df1 <- full[[condition1]]
  df2 <- full[[condition2]]
  
  df1_trunc <- na.omit(df1[match(row.names(df2), row.names(df1)), ])
  df2_trunc <- na.omit(df2[match(row.names(df1), row.names(df2)), ])
  
  if(!isTRUE(all.equal.character(rownames(df1_trunc),rownames(df2_trunc)))){
    stop("Error finding common celltypes")
  }
  
  # Declar a vector of the common celltypes
  celltypes <- row.names(df1_trunc)
  
  # Rename
  gene_df1 <- df1_trunc
  gene_df2 <- df2_trunc
  
  # Loop through celltypes
  for (item in celltypes){
    print(item)
    
    # Iterate through both directions of differential expression
    for (direction in up_down){
      
      if (direction=="UP"){
        up_degs1 <- colnames(gene_df1[item,(gene_df1[item,] > 0),drop=FALSE]) # Gets list of upregulated genes for cell type of interest for cond1
        up_degs2 <- colnames(gene_df2[item,(gene_df2[item,] > 0),drop=FALSE]) # Gets list of upregulated genes for cell type of interest for cond2
        
        delta_up <- setdiff(up_degs1,up_degs2) # Get elements in up_degs1 but not in up_degs2
        
        # Send the delta DEGs for the celltype in the direction of interest to enrichR
        if (!identical(delta_up, character(0))){
          enrich_fun(conditions,item,delta_up,direction,dbs)
        } 
      }
      
      if (direction=="DOWN"){
        down_degs1 <- colnames(gene_df1[item,(gene_df1[item,] < 0),drop=FALSE]) # Gets list of downregulated genes for cell type of interest for cond1
        down_degs2 <- colnames(gene_df2[item,(gene_df2[item,] < 0),drop=FALSE]) # Gets list of downregulated genes for cell type of interest for cond2
        
        delta_down <- setdiff(down_degs1,down_degs2) # Get elements in down_degs1 but not in down_degs2
        
        # Send the delta DEGs for the celltype in the direction of interest to enrichR
        if (!identical(delta_down, character(0))){
          enrich_fun(conditions,item,delta_down,direction,dbs)
        } 
      }
    }
  }
}

# Create object containing the DEG lists
gene_df_S <- read.csv(file =  paste0(DE_dir,"Severe_intrasev_DEG[5].csv"),check.names=FALSE,row.names=1)
gene_df_M <- read.csv(file =  paste0(DE_dir,"Mild_intrasev_DEG[5].csv"),check.names=FALSE,row.names=1)
gene_df_HC <- read.csv(file = paste0(DE_dir,"Healthy Control_intrasev_DEG[5].csv"),check.names=FALSE,row.names=1)

full <- list(gene_df_S, gene_df_M,gene_df_HC)
names(full) <- c("Severe", "Mild", "Healthy Control")

# Perform enrichR
# Make sure that the more severe condition comes first!
conditions <- c("Severe","Mild")
enrich_preprocess(conditions,full)

conditions <- c("Severe","Healthy Control")
enrich_preprocess(conditions,full)

conditions <- c("Mild","Healthy Control")
enrich_preprocess(conditions,full)

