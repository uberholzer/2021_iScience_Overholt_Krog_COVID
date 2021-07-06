# Intra-Severity Cross-Compartmental Comparison [5]
# DEGs between (celltype x, condition y) in BALF and (celltype x, condition y) in PBMC

### This notebook is used to investigate all genes found to be differentially expressed in the various patient populations when looking at scRNA-seq data of COVID patients using BALF and blood samples.

### What is needed for running:
### 1. CSV files containing the output of differential expression analysis that was carried out previously

import numpy as np
import pandas as pd
import scipy as sp
import statistics
import matplotlib
import seaborn as sns
from matplotlib import pyplot
import os
import plotly
import plotly.graph_objects as go

# Get the directories for reading/saving files
wd = os.getcwd()
print(wd)

res_dir = wd + '/Differential_Expression/DE_Plots[5]/'
DE_dir = wd + '/Differential_Expression/DE_Results[5]/'

# Function to filter genes from a DE result list    
def filter_genes(data):
    # Following DE remove genes we don't care about that could be artefactual
    # This will help to make more valid compartmental comparisons
    
    # Filter genes based on regex
    to_delete=data.filter(regex="^RP[SL]|^HB[^(P)]|^MALAT1$|^NEAT1$|^PPBP$|^MT-|^RNA18S|^RNA28S|^IGH|^IGL|^IGK|^TRAV|^TRBV|^TRDV|^TRGV", axis=0)
    
    data = data[~data.index.isin(to_delete.index)]

    return(data)

# Function getting a list of significant DEGs from a DE result list
def get_degs(condition,celltypes):
    fracs= pd.DataFrame(index=celltypes,columns = ['Frac DE'])
    degs_tot = set() # This records the TOTAL number of DEGs coming out of all celltypes

    for ctype in celltypes:
        data_unfilt = pd.read_csv(DE_dir+condition+"_"+ctype+"_DE_nCount_RNA_study.csv",index_col=0)
        data=filter_genes(data_unfilt) # Filter out ribosomal, mito, lncRNAs, and Igs              # This line debugs fine
        data = data[(data['pct.1'] > 0.1) | (data['pct.2'] > 0.1)] # Subset on min.pct       # This line debugs fine
        pass_filter = data.loc[(abs(data['avg_log2FC'])>1) & (data['p_val_adj']<.01)].index
        degs_tot.update(set(pass_filter))
        
        # Calculate fraction differentially expressed
        fraction=len(pass_filter)/data.shape[0] # or just data?
        print(condition)
        print(ctype)
        print(len(pass_filter))
        fracs.loc[ctype,'Frac DE'] = fraction
        
    return (degs_tot,fracs)


# Function saving .csv of only significant DEGs for each condition
# Can load this into GO downstream
def matches2dataframe(condition,celltypes,degs):
    combinedDF = pd.DataFrame(index = celltypes,columns=degs)
    for item in degs:
        for ctype in celltypes:
            data = pd.read_csv(DE_dir+condition+"_"+ctype+"_DE_nCount_RNA_study.csv",index_col=0)
            data=filter_genes(data) # Filter out ribosomal, mito, lncRNAs, and Igs
            data = data[(data['pct.1'] > 0.1) | (data['pct.2'] > 0.1)] # Subset on min.pct
            if item in data.index and (data.loc[item]['p_val_adj']<.01) and (abs(data.loc[item]['avg_log2FC'])>1):
                combinedDF.loc[ctype,item] = data.loc[item,'avg_log2FC']

    combinedDF = combinedDF.fillna(0)
    combinedDF.to_csv(res_dir+condition+'_intrasev_DEG[5].csv')

    
# ## Nightingale Rose Plot
# See https://medium.com/@abhishekdas.nitc/nightingale-plots-in-python-using-plotly-da42bc18d15d
# https://plotly.com/python/polar-chart/#polar-bar-chart
# https://plotly.com/python/static-image-export/

def roseplot(df,comparison):

    fig = go.Figure()

    color_P='rgb(31,158,137)'

    fig.add_trace(go.Barpolar(
        r=list(df['Frac DE']),
        theta=list(df.index),
        name='Fraction Differentially Expressed Genes',
        marker_color=color_P,
        #marker_line_color="black",
        marker_line_width=1,
        opacity=0.8))
    
    fig.update_traces(text=df.index)
    fig.update_layout(
        title=comparison+' Intra-Severity [5]',
        template=None, # plotly_dark
        font_size=16,
        legend_font_size=20,
        width=1200,
        height=800,
        polar_angularaxis_rotation=-80,

        polar = dict(
            radialaxis = dict(range=[0, 0.06], showticklabels=True, ticks='outside',linewidth = 1,gridwidth = 2),
            angularaxis = dict(showticklabels=True, ticks='outside',linewidth = 2,gridwidth = 2)
        )
    )

    fig.show()
    plotly.offline.plot(fig, filename=res_dir+comparison+'_intrasev[5].html')
    fig.write_image(res_dir+comparison+'_intrasev[5].pdf') # Must have Kaleido installed
    

#-------------------------------------------------------------------------------------------------------------------------------#    
# All DEGs Intra-Severity Compartmental Comparisons Rose Plots [5]

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
(degs,fracs) = get_degs('Severe', celltypes_S)
print(fracs)
roseplot(fracs,'Severe')
matches2dataframe('Severe',celltypes_S,degs)

# Mild
(degs,fracs) = get_degs('Mild', celltypes_M)        
print(fracs)
roseplot(fracs,'Mild')
matches2dataframe('Mild',celltypes_M,degs)

# Healthy Control
(degs,fracs) = get_degs('Healthy Control', celltypes_HC)        
print(fracs)
roseplot(fracs,'Healthy Control')
matches2dataframe('Healthy Control',celltypes_HC,degs)