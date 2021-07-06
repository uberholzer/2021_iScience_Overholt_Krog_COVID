# # GSEA for Figure 3 - Heatmaps, Dotplots, Rose Plots
# 
# ### See https://stackoverflow.com/questions/59381273/heatmap-with-circles-indicating-size-of-population for information
# ### Or see corrplot function in R

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.collections import PatchCollection
import random
import os
import plotly.graph_objects as go
import plotly

random.seed(10097)

# Get the directories for reading/saving files
wd = os.getcwd()
print(wd)


res_dir = wd + '/GSEA/GSEA_rose_dotplot/'
data_dir = wd +'/GSEA/GSEA Results/'

# ### Plotting function for Dotplot

def dotplot(pw_by_cells,p_vals,size,option):
    fig=fig = figure(num=None, figsize=size, dpi=300)
    ax  = fig.add_subplot(111)

    ylabels = pw_by_cells.index
    xlabels = pw_by_cells.columns

    # M by N matrix
    # Note M and N are backwards in the tutorial given
    
    M=len(pw_by_cells.index) # rows (pathways)
    N=len(pw_by_cells.columns) # cols (celltype)
    x, y = np.meshgrid(np.arange(N), np.arange(M))

    arr= pw_by_cells.to_numpy(dtype=float).flatten()

    s = abs(p_vals.to_numpy(dtype=float)) # This is an M by N array representing p-values
    R = s/s.max()/2.5
    
    for m in range(M):
        for n in range(N):
            if (s[m,n]<(-1*np.log10(0.05))):
                R[m,n]=0.05 # 0.5 is the "min" circle radius
        

    circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]
    col = PatchCollection(circles, array=arr, cmap="coolwarm")
    col.set_clim([-3, 3]) # This is where you set the range of the colorbar
    ax.add_collection(col)

    ax.set(xticks=np.arange(N), yticks=np.arange(M),xticklabels=xlabels, yticklabels=ylabels)
    ax.set_xticks(np.arange(N+1)-0.5, minor=True) # Need to set minor ticks to get it to work properly
    ax.set_yticks(np.arange(M+1)-0.5, minor=True)
    ax.tick_params(which='minor', length=0) # Get rid of the minor ticks
    plt.setp(ax.get_xticklabels(), rotation=90)

    #ax.grid(which='minor')

    cbar = ax.figure.colorbar(col, ax=ax)
    cbar.ax.tick_params(labelsize=8)  # set tick font size here
    cbar.set_label(label='Normalized Enrichment Score',size=14)
    
    fig.tight_layout()
    fig.savefig(res_dir+'gsea_dotplot_'+option+'_intrasev[5].pdf')

    
# ### Function for reading .csv files and choosing top GSEA pathways

def reader(filenames,celltypes,cols_rename):
    master_list=[]

    for filename in filenames:
        # Create a df for a specific cell type
        df=pd.read_csv(filename,sep='\t') 
        # Scan through the df and get the statistically significant pathways
        sigs = df[df['padj'] < 0.05] 
        temp_list = sigs['pathway'].tolist()
        # Add those to a master list
        master_list.extend(temp_list) 
#         print(master_list)
#     print('Starting num pathways: ' + str(len(master_list)))
        master_list = list(set(master_list)) # Remove duplicates from list 
        # note that the set command doesn't preserve ordering information)
    master_list=sorted(master_list,reverse=True) # alphabetize masterlist to preserve the order
#     print('Num unique pathways: '+str(len(master_list)))
#     print(filenames)
#     print(master_list)
    
    # Using the master list go back through the dataframes and extract the NE values and p-values
    pw_by_cells = pd.DataFrame(columns=celltypes,index=master_list)
    p_vals=pd.DataFrame(columns=celltypes,index=master_list)
    J_df=pd.DataFrame(columns=['Jaccard'],index=celltypes)
    unique_alpha_masterlist=[]
    unique_gamma_masterlist=[]
    #print(pw_by_cells)

    for filename, cell in zip(filenames, celltypes):

        df=pd.read_csv(filename,sep='\t')  

        # Find normalized enrichment score corresponding to each pathway, celltype
        # Find p value corresponding to each pathway, celltype
        for pw in master_list:

            if(df['pathway'].str.contains(pw).any()):
                ind=df.loc[df.isin([pw]).any(axis=1)].index.tolist()[0]
                NES=df.iloc[ind]['NES']
                padj=df.iloc[ind]['padj']

            else:
                NES=0 # Give it a zero if it's not found
                padj=1

            pw_by_cells.loc[pw][cell]=NES
            p_vals.loc[pw][cell]=-1*log10(padj) 

    # Remove the prefix "HALLMARK_" from the front of each pathway
    pw_by_cells=pw_by_cells.rename(index=lambda s: s[len('HALLMARK_'):])

    # Rename the columns for better visualization (cols_rename MUST be in the same order)
    pw_by_cells.columns = cols_rename

    return(pw_by_cells,p_vals)


# Function to prepare for dotplot, heatmap, and interferon plots
# Calls dotplot, interferons, and reader
def prep_plots(option,celltypes,size):
    filenames=[]
    # Dataset should be of the form "Lee_"
    # size is (x,y)

    for cell in celltypes:
        filenames.append(data_dir+option+cell+" fgsea.tsv")      
            
    cols_rename=celltypes

    pw_by_cells,p_vals=reader(filenames,celltypes,cols_rename)
    
    # Make the dotplot
    dotplot(pw_by_cells,p_vals,size,option)
    
    return(p_vals,pw_by_cells)


#-------------------------------------------------------------------------------------------------------------------------------#    

# Define celltypes
celltypes_HC = [
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
                "Proliferating T"]

celltypes_M = [
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
                "Proliferating T"]

celltypes_S = [
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
                "Proliferating T"]

# Severe
p_vals_BALF,pw_by_cells_BALF=prep_plots('Severe_',celltypes_S,(7,10))

# Mild
p_vals_BALF,pw_by_cells_BALF=prep_plots('Mild_',celltypes_M,(7,10))

# Healthy Control
p_vals_BALF,pw_by_cells_BALF=prep_plots('Healthy Control_',celltypes_HC,(7,10))