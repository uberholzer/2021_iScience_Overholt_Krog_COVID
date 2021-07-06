import pandas as pd
import scipy as sp
import seaborn as sns
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *
import numpy as np


# Get the directories for reading/saving files
wd = os.getcwd()

res_dir = wd + '/NicheNet/'
plots_dir = wd + '/NicheNet/Plots/'
unsub_dir = wd + '/NicheNet/Unsubtracted/'
DE_dir = wd + '/Differential_Expression/DE_Results[1]/'

def heatmap(celltypes,sample_name,sender_file_head):
    
    print(sample_name) # Progress bar
    
    # Read in the NicheNet results
    # Ligands are genes that come up in any of the cells of the receiver population for a given comparison
    ligands = set()
    for celltype in celltypes:
        data = pd.read_csv(res_dir+sample_name+celltype+'_ligand_table.txt', sep="\t",index_col=1)
        top_data = data.loc[(data['pearson']>0)].index
        ligands.update(top_data)
    ligands = list(ligands)
    df = pd.DataFrame(index = celltypes,columns = ligands)

    # Get differentially expressed ligands (DELs)
    DELs=list()
    sender_fnames = [x for x in os.listdir(DE_dir) if x.startswith(sender_file_head)] 
    for file in sender_fnames:
        data = pd.read_csv(DE_dir+file,index_col=0)
        for ligand in ligands:
            if ligand in data.index:
                if (data.loc[ligand,"p_val_adj"] < .05) and (data.loc[ligand,"avg_log2FC"]>0):
                    # Add the DEL to the list
                    DELs.append(ligand)
                    # Now print whether the DEL went up or down (use this info in the paper schematic)
                    if (data.loc[ligand,"avg_log2FC"]>0):
                        print("UP: "+ligand)
                        print(str(data.loc[ligand,"avg_log2FC"])+"\n")
                    if (data.loc[ligand,"avg_log2FC"]<0):
                        print("DOWN: "+ligand)
                        print(str(data.loc[ligand,"avg_log2FC"])+"\n")
    DELs=list(set(DELs))                

    # Fill the dataframe using the NicheNet results
    for celltype in celltypes:
        data = pd.read_csv(res_dir+sample_name+celltype+'_ligand_table.txt', sep="\t",index_col=1)
        for ligand in DELs:
            if ligand in data.loc[(data['pearson']>0)].index:
                df.loc[celltype,ligand] = data.loc[ligand,"pearson"]
    df = df.fillna(0)
    
    # Get columns where a gene has Pearson >0.1 in over 1/3 of the compartment's cell types
    df_top = df > 0.08
    df_top = df_top.sum(axis=0)
    inds_top=df_top[df_top>ceil(df.shape[0]/4)].index

    df_broad=df[inds_top]
    df_broad=df_broad.sort_index(ascending=True, axis=1) # Alphabetize columns
    
    # Print the figure
    plt.figure() # figsize=(width,height) as an argument
    sns.set(font_scale=1) 
    sns.heatmap(df_broad, annot=False,  cmap= 'PuRd', cbar_kws={'label': "Pearson Correlation"},vmin=0,vmax=0.25,square=True)
    plt.title(sample_name, fontsize =15)

    plt.tight_layout()
    plt.savefig(unsub_dir+sample_name+'BROAD.pdf',bbox='tight')
    return(df_broad)




def subtractor(df_intra,df_cross,pair_name):
    
    # Need special cases if either of these dataframes is completely empty
    if (df_intra.empty) and (not df_cross.empty):
        gray_df=0*df_cross
        df_sub=-1*df_cross
        
    if (df_cross.empty) and (not df_intra.empty):
        gray_df=0*df_intra
        df_sub=df_intra
        
    if (df_intra.empty) and (df_cross.empty):
        print("Error: Both dataframes were empty!")
        quit()
        
    if (not df_intra.empty) and (not df_cross.empty):
    # Get the columns where the subtraction will be zero
        truth_df=df_intra.combine(df_cross,np.equal,fill_value=0) # Dataframe that fills exactly common columns with 1s, 0 o/w
        gray_df=df_intra.combine(df_cross,np.maximum,fill_value=0).multiply(truth_df).T # This dataframe should be 1s when found in both compartments, 0 o/w
        
        # Subtract the dataframes
        # Fill values with zero if they don't exist
        # Columns will be the union of both dataframes
        df_sub=df_intra.combine(df_cross, pd.Series.sub, fill_value=0).T # Makes positive vals when 'intra' only, negative vals when 'cross' only
        # Negative is green meaning cross-compartmental
        # Positive is brown meaning intra-comaprtmental
    
    if not (gray_df.empty):
        plt.figure(figsize=(6,10)) # figsize=(width,height) as an argument
        sns.set(font_scale=0.5) 
        sns.heatmap(gray_df, annot=False,  cmap= 'binary', cbar_kws={'label': "Pearson's r"},square=True,vmin=0,vmax=0.4)
        plt.title(pair_name, fontsize =20)
        plt.tight_layout()
        plt.savefig(plots_dir+pair_name+'_BROAD_bw.pdf', bbox='tight')
        
#     if not (truth_df.empty):
#         plt.figure(figsize=(6,10)) # figsize=(width,height) as an argument
#         sns.set(font_scale=1) 
#         sns.heatmap(truth_df, annot=False,  cmap= 'binary', cbar_kws={'label': "Pearson's r"},square=True,vmin=0,vmax=1)
#         plt.title(pair_name, fontsize =20)
#         plt.tight_layout()
#         plt.savefig(plots_dir+pair_name+'_BROAD_TRUTH_bw.pdf', bbox='tight')
    
    # Print the figure
    if not (df_sub.empty):
        plt.figure(figsize=(6,10)) # figsize=(width,height) as an argument
        sns.set(font_scale=0.5) 
        sns.heatmap(df_sub, annot=False,  cmap= 'BrBG_r', cbar_kws={'label': "Pearson's r"},square=True,center=0,vmin=-0.2,vmax=0.2)
        plt.title(pair_name, fontsize =20)

        plt.tight_layout()
        plt.savefig(plots_dir+pair_name+'_BROAD_color.pdf', bbox='tight')




# Define the cell types for each compartmental comparison in the right order

receivers_liao_SC=["AM",
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
                  ]

receivers_liao_SM=["AM",
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
                    "Proliferating T"]

receivers_liao_MC=["AM",
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
                  ]

receivers_lee_SC=["CD14 Mono",
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
                   "Proliferating T"]

receivers_lee_SM=["CD14 Mono",
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
#                    "Plasma cell",
                   "Proliferating T"]

receivers_lee_MC=receivers_lee_SM




# Loop through the various compartmental comparisons
# Generate heatmap for each (receiver, comparison) pair

senders=["Liao","Lee"]
receivers=["Liao","Lee"]
comparisons=["SM","SC","MC"]

for receiver in receivers:
    for comparison in comparisons:
        pair_name=receiver+"_"+comparison+"" 
        for sender in senders:
            if comparison=="SM":
                comp_txt="_Severe_Mild_"
                if receiver=="Liao":
                    celltypes=receivers_liao_SM
                if receiver=="Lee":
                    celltypes=receivers_lee_SM
            if comparison=="SC":
                comp_txt="_Severe_Healthy Control_"
                if receiver=="Liao":
                    celltypes=receivers_liao_SC
                if receiver=="Lee":
                    celltypes=receivers_lee_SC
            if comparison=="MC":
                comp_txt="_Mild_Healthy Control_"
                if receiver=="Liao":
                    celltypes=receivers_liao_MC
                if receiver=="Lee":
                    celltypes=receivers_lee_MC
            
                    
            sample_name=sender+"_to_"+receiver+"_"+comparison+"_"
            
            sender_file=sender+comp_txt       
            
            # generate heatmaps 1 and 2
            if sender==receiver:
                df_intra=heatmap(celltypes,sample_name,sender_file)
            if sender!=receiver:
                df_cross=heatmap(celltypes,sample_name,sender_file)
            
        # combine heatmaps for a given (receiver, comparison) pair
        subtractor(df_intra,df_cross,pair_name)