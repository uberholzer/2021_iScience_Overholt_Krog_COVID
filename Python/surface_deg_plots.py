# Surface Protein DEG Investigation

# To convert from .ipynb notebook to .py, run this in the notebook
# get_ipython().system('jupyter nbconvert --to script surface_proteins.ipynb')

### This notebook is used to investigate the surface proteins that are found to be differentially expressed in the various patient populations when looking at scRNA-seq data of COVID patients using BALF and blood samples.

### What is needed for running:
### 1. Excel speadsheet of the surface protein atlas downloaded from https://wlab.ethz.ch/cspa/
### 2. CSV files containing the output of differential expression analysis that was carried out previously

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

res_dir = wd + '/Differential_Expression/DE_Plots[1_6]/Surface_markers/'
DE_dir = wd + '/Differential_Expression/DE_Results[1]/'

# First we import the database of known surface proteins and take only those labeled as "high confidence"
surface_prot = pd.read_excel("human_surfaceproteins.xlsx",index_col = 7)
surf_prot_high = surface_prot.loc[surface_prot['CSPA category']=='1 - high confidence'].index


# ## Functions for finding surface proteins and plotting

# Get DE proteins that match the surface protein atlas
def get_degs(prefix,condition,celltypes):
    matches = set()
    for ctype in celltypes:
        data = pd.read_csv(DE_dir+prefix+condition+ctype+"_DE_nCount_RNA.csv",index_col=0) # Read file
        data=filter_genes(data) # Filter out ribosomal, mito, lncRNAs, and Igs
        data = data[(data['pct.1'] > 0.1) | (data['pct.2'] > 0.1)] # Subset here based on min.pct
        pass_filter = data.loc[(abs(data['avg_log2FC'])>1) & (data['p_val_adj']<.01)].index
        match = set(surf_prot_high) & set(pass_filter)
        matches.update(match)
    return matches

# Function saving .csv of only significant DEGs
def matches2dataframe(prefix,condition,celltypes,ordered_celltypes,matches):
    combinedDF = pd.DataFrame(index = ordered_celltypes,columns=matches)
    for item in matches:
        for ctype in celltypes:
            data = pd.read_csv(DE_dir+prefix+condition+ctype+"_DE_nCount_RNA.csv",index_col=0) # Read file
            data=filter_genes(data) # Filter out ribosomal, mito, lncRNAs, and Igs
            data = data[(data['pct.1'] > 0.1) | (data['pct.2'] > 0.1)] # Subset here based on min.pct
            if item in data.index and (data.loc[item]['p_val_adj']<.01) and (abs(data.loc[item]['avg_log2FC'])>1):
                combinedDF.loc[ctype,item] = data.loc[item,'avg_log2FC']

            
    combinedDF = combinedDF.fillna(0)
    return combinedDF


def plotMatches(dataframe,comparison,sample_type):
    fig = pyplot.figure(figsize=(45, 20))# width and height in inches
    sns.set(font_scale=1.55) 
    g = sns.heatmap(dataframe, annot=False,  cmap= 'RdBu_r', center=0, cbar_kws={'label': "Log$_2$ Fold Change"})
    g.figure.axes[-1].yaxis.label.set_size(30)
    g.set_yticklabels(g.get_ymajorticklabels(), fontsize = 25)
    pyplot.title(comparison + ' Surface Markers Differential Expression in '+sample_type,fontsize= 30,weight='bold')
    fig.tight_layout()
    pyplot.savefig(res_dir+sample_type + comparison+'SURFACE_MARKERS.pdf', dpi=300)
    dataframe.to_csv(res_dir+sample_type + comparison+'SURFACE_MARKERS.csv')
    return


# ## Functions for making compartmental comparisons
def jaccard(BALF_dataset,PBMC_dataset,comparison):
    
    pbmc_markers = pd.read_csv(res_dir+PBMC_dataset+'_'+comparison+'_SURFACE_MARKERS.csv',index_col=0)
    balf_markers = pd.read_csv(res_dir+BALF_dataset+'_'+comparison+'_SURFACE_MARKERS.csv',index_col=0)
    
    overlap_celltypes=list(set(balf_markers.index) & set(pbmc_markers.index))
    print(overlap_celltypes)

    df = pd.DataFrame(columns = overlap_celltypes,index=['Frac BALF Only','Frac PBMC Only','Jaccard Contradirectional','Jaccard Codirectional'])

    for ctype in overlap_celltypes:
        
        # Filter for only nonzero columns in each BALF cell type, count how many
        # Create separate lists for upreg and downreg markers
        foo=balf_markers.loc[ctype]
        balf_list_up=foo[(foo!=0)&(foo>0)].index
        balf_list_down=foo[(foo!=0)&(foo<0)].index
        n_balf_up=len(balf_list_up)
        n_balf_down=len(balf_list_down)
        n_balf=n_balf_up+n_balf_down

        # Filter for only nonzero columns in each PBMC cell type, count how many
        # Create separate lists for upreg and downreg markers
        foo=pbmc_markers.loc[ctype]
        pbmc_list_up=foo[(foo!=0)&(foo>0)].index
        pbmc_list_down=foo[(foo!=0)&(foo<0)].index
        n_pbmc_up=len(pbmc_list_up)
        n_pbmc_down=len(pbmc_list_down)
        n_pbmc=n_pbmc_up+n_pbmc_down

        # Calculate intersection
        intersect_genes_co=list(set(balf_list_up) & set(pbmc_list_up)) + list(set(balf_list_down) & set(pbmc_list_down))
        intersect_genes_contra=list(set(balf_list_up) & set(pbmc_list_down)) + list(set(balf_list_down) & set(pbmc_list_up))
        
        U=n_pbmc+n_balf-len(intersect_genes_co)-len(intersect_genes_contra)
        I_co=len(intersect_genes_co)
        I_contra=len(intersect_genes_contra)
        
        # Save BALF and PBMC fraction
        # Calculate Jaccard index
        
        if (U!=0):
            df.loc['Frac BALF Only',ctype]=(n_balf-I_co-I_contra)/U
            df.loc['Frac PBMC Only',ctype]=(n_pbmc-I_co-I_contra)/U
            
            J_co=I_co/U
            df.loc['Jaccard Codirectional',ctype]=J_co
            
            J_contra=I_contra/U
            df.loc['Jaccard Contradirectional',ctype]=J_contra
            
        else:
            df.loc['Frac BALF Only',ctype]=0
            df.loc['Frac PBMC Only',ctype]=0
            df.loc['Jaccard Codirectional',ctype]=0
            df.loc['Jaccard Contradirectional',ctype]=0

        # Print information
        print(comparison)
        print(BALF_dataset+PBMC_dataset)
        print(ctype)
        print('n BALF: ' + str(n_balf))
        print('n PBMC: ' + str(n_pbmc))
        print('Union: '+str(U))
        print('Contra-Intersection: '+str(I_contra))
        print('Co-Intersection: '+str(I_co))
        if (I_co!=0):
            print(intersect_genes_co) # Show the intersecting genes if any
        print("\n")

    df=df.T
    return(df)

# Function for making within compartment comparisons (control [6])
def jaccard_ctrl(dataset1,dataset2,comparison):
    
    ds2_markers = pd.read_csv(res_dir+dataset2+'_'+comparison+'_SURFACE_MARKERS.csv',index_col=0)
    ds1_markers = pd.read_csv(res_dir+dataset1+'_'+comparison+'_SURFACE_MARKERS.csv',index_col=0)
    
    overlap_celltypes=list(set(ds1_markers.index) & set(ds2_markers.index))
    print(overlap_celltypes)

    df = pd.DataFrame(columns = overlap_celltypes,index=['Frac ds1 Only','Frac ds2 Only','Jaccard Contradirectional','Jaccard Codirectional'])

    for ctype in overlap_celltypes:
        
        # Filter for only nonzero columns in each ds1 cell type, count how many
        # Create separate lists for upreg and downreg markers
        foo=ds1_markers.loc[ctype]
        ds1_list_up=foo[(foo!=0)&(foo>0)].index
        ds1_list_down=foo[(foo!=0)&(foo<0)].index
        n_ds1_up=len(ds1_list_up)
        n_ds1_down=len(ds1_list_down)
        n_ds1=n_ds1_up+n_ds1_down

        # Filter for only nonzero columns in each ds2 cell type, count how many
        # Create separate lists for upreg and downreg markers
        foo=ds2_markers.loc[ctype]
        ds2_list_up=foo[(foo!=0)&(foo>0)].index
        ds2_list_down=foo[(foo!=0)&(foo<0)].index
        n_ds2_up=len(ds2_list_up)
        n_ds2_down=len(ds2_list_down)
        n_ds2=n_ds2_up+n_ds2_down

        # Calculate intersection
        intersect_genes_co=list(set(ds1_list_up) & set(ds2_list_up)) + list(set(ds1_list_down) & set(ds2_list_down))
        intersect_genes_contra=list(set(ds1_list_up) & set(ds2_list_down)) + list(set(ds1_list_down) & set(ds2_list_up))
        
        U=n_ds2+n_ds1-len(intersect_genes_co)-len(intersect_genes_contra)
        I_co=len(intersect_genes_co)
        I_contra=len(intersect_genes_contra)
        
        # Save ds1 and ds2 fraction
        # Calculate Jaccard index
        
        if (U!=0):
            df.loc['Frac Dataset1 Only',ctype]=(n_ds1-I_co-I_contra)/U
            df.loc['Frac Dataset2 Only',ctype]=(n_ds2-I_co-I_contra)/U
            
            J_co=I_co/U
            df.loc['Jaccard Codirectional',ctype]=J_co
            
            J_contra=I_contra/U
            df.loc['Jaccard Contradirectional',ctype]=J_contra
            
        else:
            df.loc['Frac ds1 Only',ctype]=0
            df.loc['Frac ds2 Only',ctype]=0
            df.loc['Jaccard Codirectional',ctype]=0
            df.loc['Jaccard Contradirectional',ctype]=0

        # Print information
        print(ctype)
        print('n ds1: ' + str(n_ds1))
        print('n ds2: ' + str(n_ds2))
        print('Union: '+str(U))
        print('Contra-Intersection: '+str(I_contra))
        print('Co-Intersection: '+str(I_co))
        if (I_co!=0):
            print(intersect_genes_co) # Show the intersecting genes if any
        print("\n")

    df=df.T
    return(df)


# ## Function for creating Nightingale Rose Plots

# See https://medium.com/@abhishekdas.nitc/nightingale-plots-in-python-using-plotly-da42bc18d15d
# https://plotly.com/python/polar-chart/#polar-bar-chart
# https://plotly.com/python/static-image-export/

def roseplot(df,sets,comparison):

    fig = go.Figure()

    color_J_co='rgb(253,231,37)'
    color_J_contra='rgb(109,205,89)'
    color_P='rgb(31,158,137)'
    color_B='rgb(68,1,84)'

    fig.add_trace(go.Barpolar(
        r=list(df['Frac BALF Only']),
        theta=list(df.index),
        name='BALF Only',
        marker_color=color_B,
        #marker_line_color="black",
        marker_line_width=1,
        opacity=0.8
    ))
    fig.add_trace(go.Barpolar(
        r=list(df['Frac PBMC Only']),
        theta=list(df.index),
        name='PBMC Only',
        marker_color=color_P,
        #marker_line_color="black",
        marker_line_width=1,
        opacity=0.8
    ))
    fig.add_trace(go.Barpolar(
        r=list(df['Jaccard Contradirectional']),
        theta=list(df.index),
        name='Jaccard Contradirectional',
        marker_color=color_J_contra,
        #marker_line_color="black",
        marker_line_width=1,
        opacity=0.8
    ))
    fig.add_trace(go.Barpolar(
        r=list(df['Jaccard Codirectional']),
        theta=list(df.index),
        name='Jaccard Codirectional',
        marker_color=color_J_co,
        #marker_line_color="black",
        marker_line_width=1,
        opacity=0.8
    ))


    # fig.update_traces(text=df.index)
    fig.update_layout(
        title=sets+' '+comparison+' Surface Markers',
        template=None, # plotly_dark
        font_size=16,
        legend_font_size=20,
        width=1200,
        height=800,
        polar_angularaxis_rotation=-80,

        polar = dict(
            radialaxis = dict(range=[0, 1.2], showticklabels=False, ticks='outside',linewidth = 1,gridwidth = 2),
            angularaxis = dict(showticklabels=True, ticks='outside',linewidth = 2,gridwidth = 2)
        )
    )


    fig.show()
    plotly.offline.plot(fig, filename=res_dir+'rose[1]_SURFACE_MARKERS_'+sets+'_'+comparison+'.html')
    fig.write_image(res_dir+'rose[1]_SURFACE_MARKERS_'+sets+'_'+comparison+'.pdf') # Must have Kaleido installed
    
def roseplot_ctrl(df,sets,comparison):

    fig = go.Figure()

    color_J_co='rgb(253,231,37)'
    color_J_contra='rgb(109,205,89)'
    color_P='rgb(31,158,137)'
    color_B='rgb(68,1,84)'

    fig.add_trace(go.Barpolar(
        r=list(df['Frac Dataset1 Only']),
        theta=list(df.index),
        name='Dataset1 Only',
        marker_color=color_B,
        #marker_line_color="black",
        marker_line_width=1,
        opacity=0.8
    ))
    fig.add_trace(go.Barpolar(
        r=list(df['Frac Dataset2 Only']),
        theta=list(df.index),
        name='Dataset2 Only',
        marker_color=color_P,
        #marker_line_color="black",
        marker_line_width=1,
        opacity=0.8
    ))
    fig.add_trace(go.Barpolar(
        r=list(df['Jaccard Contradirectional']),
        theta=list(df.index),
        name='Jaccard Contradirectional',
        marker_color=color_J_contra,
        #marker_line_color="black",
        marker_line_width=1,
        opacity=0.8
    ))
    fig.add_trace(go.Barpolar(
        r=list(df['Jaccard Codirectional']),
        theta=list(df.index),
        name='Jaccard Codirectional',
        marker_color=color_J_co,
        #marker_line_color="black",
        marker_line_width=1,
        opacity=0.8
    ))


    # fig.update_traces(text=df.index)
    fig.update_layout(
        title=sets+' '+comparison+' Surface Markers',
        template=None, # plotly_dark
        font_size=16,
        legend_font_size=20,
        width=1200,
        height=800,
        polar_angularaxis_rotation=-80,

        polar = dict(
            radialaxis = dict(range=[0, 1.2], showticklabels=False, ticks='outside',linewidth = 1,gridwidth = 2),
            angularaxis = dict(showticklabels=True, ticks='outside',linewidth = 2,gridwidth = 2)
        )
    )

    
    fig.show()
    plotly.offline.plot(fig, filename=res_dir+'rose[6]_SURFACE_MARKERS_'+sets+'_'+comparison+'.html')
    fig.write_image(res_dir+'rose[6]_SURFACE_MARKERS_'+sets+'_'+comparison+'.pdf') # Must have Kaleido installed
    
    
def filter_genes(data):
    # Following DE remove genes we don't care about that could be artefactual
    # This will help to make more valid compartmental comparisons
    
    # Filter genes based on regex
    to_delete=data.filter(regex="^RP[SL]|^HB[^(P)]|^MALAT1$|^NEAT1$|^PPBP$|^MT-|^RNA18S|^RNA28S|^IGH|^IGL|^IGK|^TRAV|^TRBV|^TRDV|^TRGV", axis=0)
    
    data = data[~data.index.isin(to_delete.index)]

    return(data)
#-------------------------------------------------------------------------------------------------------------------------------#    
# ## Heatmaps of surface markers

########################## LIAO BALF ##########################
dataset="Liao_BALF_"

# Severe vs. Healthy Control
celltypes = [
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

x = get_degs(prefix = "Liao_", condition = 'Severe_Healthy Control_', celltypes=celltypes)   
y = matches2dataframe(prefix = "Liao_", condition = 'Severe_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Severe_vs_Healthy_Control_",sample_type = dataset)

# Mild vs. Healthy Control (use same celltypes)
x = get_degs(prefix = "Liao_", condition = 'Mild_Healthy Control_', celltypes=celltypes)   
y = matches2dataframe(prefix = "Liao_", condition = 'Mild_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Mild_vs_Healthy_Control_",sample_type = dataset)

# Severe vs. Mild 
celltypes = [
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
              # "Plasma cell",
              "Proliferating T"]

x = get_degs(prefix = "Liao_",condition = 'Severe_Mild_', celltypes=celltypes)     
y = matches2dataframe(prefix = "Liao_", condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Severe_vs_Mild_", sample_type = dataset)


########################## WAUTERS BALF ##########################

dataset="Wauters_BALF_"

# Severe vs. Mild
celltypes = [
            "CD14 Mono",
            # "CD16 Mono",
            # "CD4 Naive T",
            "CD4 Treg",
            "CD8 Effector T",
            # "CD8 Naive T",
            # "cDC1",
            "cDC2",
            "Inflammatory MP",
            # "Intermediate Mono",
            # "Mature B",
            "Neutrophil",
            "NK",
            # "pDC",
            # "Plasma cell",
            "Proliferating T"]

x = get_degs(prefix = "Wauters_", condition = 'Severe_Mild_', celltypes=celltypes)   
y = matches2dataframe(prefix = "Wauters_", condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Severe_vs_Mild_",sample_type = dataset)

########################## WILK PBMC ##########################

dataset="Wilk_PBMC_"

# Combined Severe vs. Healthy Control
celltypes = [
              "CD14 Mono",
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
              "Proliferating T"]

x = get_degs(prefix = "Wilk_",condition = 'Severe_Healthy Control_', celltypes=celltypes)    
y = matches2dataframe(prefix = "Wilk_", condition = 'Severe_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Severe_vs_Healthy_Control_",sample_type = dataset)


########################## LEE PBMC ##########################

dataset="Lee_PBMC_"
prefix="Lee_"

# Severe vs. Healthy Control
celltypes = [
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
                "Proliferating T"]

x = get_degs(prefix = prefix,condition = 'Severe_Healthy Control_', celltypes=celltypes)    
y = matches2dataframe(prefix = prefix, condition = 'Severe_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Severe_vs_Healthy_Control_",sample_type = dataset)

# Severe vs. Mild
x = get_degs(prefix = prefix,condition = 'Severe_Mild_', celltypes=celltypes)               
y = matches2dataframe(prefix = prefix, condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Severe_vs_Mild_",sample_type = dataset)

# Mild vs. Healthy Control
x = get_degs(prefix = prefix,condition = 'Mild_Healthy Control_', celltypes=celltypes)               
y = matches2dataframe(prefix = prefix, condition = 'Mild_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Mild_vs_Healthy_Control_",sample_type = dataset)


########################## ARUN PBMC ##########################

dataset="Arun_PBMC_"
prefix="Arun_"

# Severe vs. Healthy Control
celltypes = [
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
                "Plasma cell",
                "Proliferating T"]

x = get_degs(prefix = prefix,condition = 'Severe_Healthy Control_', celltypes=celltypes)             
y = matches2dataframe(prefix = prefix, condition = 'Severe_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Severe_vs_Healthy_Control_",sample_type = dataset)

# Severe vs. Mild
x = get_degs(prefix = prefix,condition = 'Severe_Mild_', celltypes=celltypes)        
y = matches2dataframe(prefix = prefix, condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Severe_vs_Mild_",sample_type = dataset)

# Mild vs. Healthy Control
x = get_degs(prefix = prefix,condition = 'Mild_Healthy Control_', celltypes=celltypes)        
y = matches2dataframe(prefix = prefix, condition = 'Mild_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Mild_vs_Healthy_Control_",sample_type = dataset)


########################## SCHULTE PBMC ##########################

dataset="Schulte_PBMC_"
prefix="Schulte_"

# Severe vs. Healthy Control
celltypes = [
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
                "Proliferating T"]                


# Severe vs. Mild
x = get_degs(prefix = prefix,condition = 'Severe_Mild_', celltypes=celltypes)        
y = matches2dataframe(prefix = prefix, condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),matches = x)
z = plotMatches(dataframe = y,comparison = "Severe_vs_Mild_",sample_type = dataset)

    
#-------------------------------------------------------------------------------------------------------------------------------#    
# # Surface Markers Compartmental Comparison Rose Plots [1]

########################## LIAO VS WILK ##########################
BALF_dataset="Liao_BALF"
PBMC_dataset="Wilk_PBMC"

# (1) Severe vs Healthy
comparison="Severe_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Liao_vs_Wilk",comparison)


########################## LIAO VS LEE ##########################
BALF_dataset="Liao_BALF"
PBMC_dataset="Lee_PBMC"

# (2) Severe vs Healthy
comparison="Severe_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Liao_vs_Lee",comparison)

# (3) Severe vs. Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Liao_vs_Lee",comparison)

# (4) Mild vs. Healthy Control
comparison="Mild_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Liao_vs_Lee",comparison)


########################## LIAO VS ARUN ##########################
BALF_dataset="Liao_BALF"
PBMC_dataset="Arun_PBMC"

# (5) Severe vs Healthy
comparison="Severe_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Liao_vs_Arun",comparison)

# (6) Severe vs. Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Liao_vs_Arun",comparison)

# (7) Mild vs. Healthy Control
comparison="Mild_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Liao_vs_Arun",comparison)


########################## LIAO VS SCHULTE ##########################
BALF_dataset="Liao_BALF"
PBMC_dataset="Schulte_PBMC"

# (8) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Liao_vs_Schulte",comparison)

########################## WAUTERS VS LEE ##########################
BALF_dataset="Wauters_BALF"
PBMC_dataset="Lee_PBMC"

# (9) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Wauters_vs_Lee",comparison)

########################## WAUTERS VS ARUN ##########################
BALF_dataset="Wauters_BALF"
PBMC_dataset="Arun_PBMC"

# (10) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Wauters_vs_Arun",comparison)

########################## WAUTERS VS SCHULTE ##########################
BALF_dataset="Wauters_BALF"
PBMC_dataset="Schulte_PBMC"

# (11) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)
roseplot(df,"Wauters_vs_Schulte",comparison)

#-------------------------------------------------------------------------------------------------------------------------------#    

# # Surface Markers Intra-Compartmental Control Rose Plots [6]

########################## LIAO VS WAUTERS ##########################
BALF_dataset1="Liao_BALF"
BALF_dataset2="Wauters_BALF"

# (1) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard_ctrl(BALF_dataset1,BALF_dataset2,comparison)
roseplot_ctrl(df,"Liao_vs_Wauters_[6]",comparison)

########################## WILK VS LEE ##########################
PBMC_dataset1="Wilk_PBMC"
PBMC_dataset2="Lee_PBMC"

# (2) Severe vs Healthy
comparison="Severe_vs_Healthy_Control"
df=jaccard_ctrl(PBMC_dataset1,PBMC_dataset2,comparison)
roseplot_ctrl(df,"Wilk_vs_Lee_[6]",comparison)

########################## WILK VS ARUN ##########################
PBMC_dataset1="Wilk_PBMC"
PBMC_dataset2="Arun_PBMC"

# (3) Severe vs Healthy
comparison="Severe_vs_Healthy_Control"
df=jaccard_ctrl(PBMC_dataset1,PBMC_dataset2,comparison)
roseplot_ctrl(df,"Wilk_vs_Arun_[6]",comparison)

########################## LEE VS ARUN ##########################
PBMC_dataset1="Lee_PBMC"
PBMC_dataset2="Arun_PBMC"

# (4) Severe vs Healthy
comparison="Severe_vs_Healthy_Control"
df=jaccard_ctrl(PBMC_dataset1,PBMC_dataset2,comparison)
roseplot_ctrl(df,"Lee_vs_Arun_[6]",comparison)

# (5) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard_ctrl(PBMC_dataset1,PBMC_dataset2,comparison)
roseplot_ctrl(df,"Lee_vs_Arun_[6]",comparison)

# (6) Mild vs Healthy
comparison="Mild_vs_Healthy_Control"
df=jaccard_ctrl(PBMC_dataset1,PBMC_dataset2,comparison)
roseplot_ctrl(df,"Lee_vs_Arun_[6]",comparison)

########################## LEE VS SCHULTE ##########################
PBMC_dataset1="Lee_PBMC"
PBMC_dataset2="Schulte_PBMC"

# (7) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard_ctrl(PBMC_dataset1,PBMC_dataset2,comparison)
roseplot_ctrl(df,"Lee_vs_Schulte_[6]",comparison)

########################## ARUN VS SCHULTE ##########################
PBMC_dataset1="Arun_PBMC"
PBMC_dataset2="Schulte_PBMC"

# (8) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard_ctrl(PBMC_dataset1,PBMC_dataset2,comparison)
roseplot_ctrl(df,"Arun_vs_Schulte_[6]",comparison)
