library(stringr)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
library(nichenetr)
library(tidyverse)
library(future)
set.seed(10097)


wd <- getwd()
res_dir<-paste0(wd,"/NicheNet")


ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction") # Use only "bona fide" high-confidence ligands & receptors
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

get_celltypes <- function(comparison,sender_set,receiver_set){

  # Define celltype lists
  celltype_liao_SC <- c("AM",
                        "Epithelial",
                        "CD14 Mono",
                        "CD16 Mono",
                        "CD4 Naive T",
                        "CD4 Treg",
                        "CD8 Effector T",
                        # "CD8 Naive T",
                        "cDC1",
                        "cDC2",
                        "Inflammatory MP",
                        # "Intermediate Mono",
                        "Mature B",
                        "Neutrophil",
                        "NK"
                        # "pDC",
                        # "Plasma cell",
                        # "Proliferating T"
                        )

  celltype_liao_SM <- c( 
                        "AM",
                        "Epithelial",
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
                        "pDC",
                        # "Plasma cell",
                        "Proliferating T"
                        )
  
  celltype_liao_MC <- c("AM",
                        "Epithelial",
                        "CD14 Mono",
                        "CD16 Mono",
                        "CD4 Naive T",
                        "CD4 Treg",
                        "CD8 Effector T",
                        # "CD8 Naive T",
                        # "cDC1",
                        "cDC2",
                        "Inflammatory MP",
                        # "Intermediate Mono",
                        "Mature B",
                        "Neutrophil",
                        "NK"
                        # "pDC",
                        # "Plasma cell",
                        # "Proliferating T"
                        )

  celltype_lee_SC <- c("CD14 Mono",
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
                       # "Plasma cell",
                       "Proliferating T") 

  celltype_lee_SM <-  c(
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
                         # "Plasma cell",
                         "Proliferating T") 
  
  celltype_lee_MC <- celltype_lee_SC


  # Choose the celltype lists for sender and receiver
  if (sender_set=="Liao" & comparison=="SM"){
    ctypes_sender <- celltype_liao_SM
  }

  if (sender_set=="Lee" & comparison=="SM"){
    ctypes_sender <- celltype_lee_SM
  }

  if (receiver_set=="Liao"& comparison=="SM"){
    ctypes_receiver <- celltype_liao_SM
  }

  if (receiver_set=="Lee" & comparison=="SM"){
    ctypes_receiver <- celltype_lee_SM
  }

  if (sender_set=="Liao" & comparison=="SC"){
    ctypes_sender <- celltype_liao_SC
  }

  if (sender_set=="Lee" & comparison=="SC"){
    ctypes_sender <- celltype_lee_SC
  }

  if (receiver_set=="Liao"& comparison=="SC"){
    ctypes_receiver <- celltype_liao_SC
  }

  if (receiver_set=="Lee" & comparison=="SC"){
    ctypes_receiver <- celltype_lee_SC
  }
  
  if (sender_set=="Liao" & comparison=="MC"){
    ctypes_sender <- celltype_liao_MC
  }
  
  if (sender_set=="Lee" & comparison=="MC"){
    ctypes_sender <- celltype_lee_MC
  }
  
  if (receiver_set=="Liao"& comparison=="MC"){
    ctypes_receiver <- celltype_liao_MC
  }
  
  if (receiver_set=="Lee" & comparison=="MC"){
    ctypes_receiver <- celltype_lee_MC
  }
  

  return(list("ctypes_sender"=ctypes_sender,"ctypes_receiver"=ctypes_receiver))
}

niche_net <- function(ctypes_sender,full_obj_sender,ctypes_receiver,full_obj_receiver,comparison,sender_set,receiver_set){

  print(paste0("Current comparison: ",comparison))
  print(paste0("Sender compartment: ",sender_set))
  print(paste0("Receiver compartment: ",receiver_set))

  # Change wd for printing
  setwd(res_dir)

  # Get conditions for DE
  # IMPORTANT: condition_inducer must be the MORE SEVERE of the two conditions (e.g. "mild" for "MC")
  if (comparison=="SM"){
    condition_inducer = "Severe"
    condition_reference = "Mild"
  }

  if (comparison=="SC"){
    condition_inducer = "Severe"
    condition_reference = "Healthy Control"
  }

  if (comparison=="MC"){
    condition_inducer = "Mild"
    condition_reference = "Healthy Control"
  }
  
  # Get the subset of sender cells from the MORE SEVERE of the conditions (the inducer subset)
  Idents(object=full_obj_sender) <- "group"
  inducerSubset_sender <- subset(full_obj_sender,idents=condition_inducer)
  
  # Set default idents and assays
  Idents(object=full_obj_receiver) <- "celltype_fine"
  DefaultAssay(full_obj_receiver) <- "RNA"
  Idents(object=inducerSubset_sender) <- "celltype_fine"
  DefaultAssay(inducerSubset_sender) <- "RNA"


  # Run the loop
  for (receiver in ctypes_receiver){

    print(paste0("Reciever celltype: ",receiver))

    expressed_genes_receiver = get_expressed_genes(receiver, full_obj_receiver, pct = 0.10,assay_oi = 'RNA')
    background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

    list_expressed_genes_sender = ctypes_sender %>% unique() %>% lapply(get_expressed_genes, inducerSubset_sender, 0.10,assay_oi = 'RNA')
    expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

    seurat_obj_receiver= subset(full_obj_receiver, idents = receiver)
    seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]])

    # Get DEGs in the receiver cell type between the populations of interest
    DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_inducer, ident.2 = condition_reference, min.pct = 0.10) %>%
      rownames_to_column("gene")

    geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.36) %>% pull(gene)
    geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

    ligands = lr_network %>% pull(from) %>% unique()
    receptors = lr_network %>% pull(to) %>% unique()

    expressed_ligands = intersect(ligands,expressed_genes_sender)
    expressed_receptors = intersect(receptors,expressed_genes_receiver)

    potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

    ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

    ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
    ligand_activities
    write.table(ligand_activities, paste0(sender_set,"_to_",receiver_set,"_",comparison,"_",receiver,"_ligand_table.txt"), append = FALSE, sep = "\t", dec = ".",
                row.names = TRUE, col.names = TRUE)
  }
}

## Return to main and run NicheNet
main_obj <- readRDS("main_fine_labeled.rds")

# Set up multicore parallelization
plan("multicore") 
plan()
options(future.globals.maxSize = 1000 * 1024^9,seed=TRUE)

# comparisons <- c("SC","SM","MC")
comparisons <- c("SM")
datasets <- c("Liao","Lee")

for (comparison in comparisons){
  for (dataset1 in datasets){
    sender_set <- dataset1

    for (dataset2 in datasets){
      receiver_set <- dataset2

      # Get celltypes for the appropriate compartmental comparison
      list <- get_celltypes(comparison,sender_set,receiver_set)

      # Extract elements from the list
      ctypes_sender=list$ctypes_sender
      ctypes_receiver=list$ctypes_receiver
      
      # Subset out the donors of interest into their own Seurat objects
      Idents(object=main_obj) <- "study"
      full_obj_sender <- subset(main_obj, idents = sender_set)
      full_obj_receiver <- subset(main_obj, idents = receiver_set)

      # Run NicheNet
      niche_net(ctypes_sender,full_obj_sender,ctypes_receiver,full_obj_receiver,comparison,sender_set,receiver_set)
    }
  }
}

print("Script finished!")