# # Differentially Expressed Genes Investigation

# ### This notebook is used to investigate all genes found to be differentially expressed in the various patient populations when looking at scRNA-seq data of COVID patients using BALF and blood samples.
# 
# ### What is needed for running:
# ### 1. CSV files containing the output of differential expression analysis that was carried out previously

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

res_dir = wd + '/Differential_Expression/DE_Plots[3]/'
DE_dir = wd + '/Differential_Expression/DE_Results[3]/'


# Function getting a list of significant DEGs
def gets_degs(prefix,condition,celltypes):
    degs = set()
    for ctype in celltypes:
        data = pd.read_csv(DE_dir+prefix+condition+ctype+"_DE_nCount_RNA_donor.csv",index_col=0)
        data=filter_genes(data) # Filter out ribosomal, mito, lncRNAs, and Igs
        data = data[(data['pct.1'] > 0.25) | (data['pct.2'] > 0.25)] # Subset on min.pct
        pass_filter = data.loc[(abs(data['avg_log2FC'])>1.5) & (data['p_val_adj']<.01)].index
        degs.update(set(pass_filter))
    return degs


# Function saving .csv of only significant DEGs
def matches2dataframe(prefix,condition,celltypes,ordered_celltypes,comparison,sample_type,degs):
    combinedDF = pd.DataFrame(index = ordered_celltypes,columns=degs)
    for item in degs:
        for ctype in celltypes:
            data = pd.read_csv(DE_dir+prefix+condition+ctype+"_DE_nCount_RNA_donor.csv",index_col=0)
            data=filter_genes(data) # Filter out ribosomal, mito, lncRNAs, and Igs
            data = data[(data['pct.1'] > 0.25) | (data['pct.2'] > 0.25)] # Subset on min.pct
            if item in data.index and (data.loc[item]['p_val_adj']<.01) and (abs(data.loc[item]['avg_log2FC'])>1.5):
                combinedDF.loc[ctype,item] = data.loc[item,'avg_log2FC']

    combinedDF = combinedDF.fillna(0)
    
    combinedDF.to_csv(res_dir+sample_type + comparison+'STRINGENT_DEG.csv')


# ## Function for making compartmental comparisons
def jaccard(BALF_dataset,PBMC_dataset,comparison):

    pbmc_markers = pd.read_csv(res_dir+PBMC_dataset+'_'+comparison+'_STRINGENT_DEG.csv',index_col=0)
    balf_markers = pd.read_csv(res_dir+BALF_dataset+'_'+comparison+'_STRINGENT_DEG.csv',index_col=0)

    overlap_celltypes=list(set(balf_markers.index) & set(pbmc_markers.index))

    df = pd.DataFrame(columns = overlap_celltypes,index=['Frac BALF Only','Frac PBMC Only','Jaccard Contradirectional','Jaccard Codirectional'])
    
    int_list_up=[]
    int_list_down=[]

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
        intersect_genes_co_up=list(set(balf_list_up) & set(pbmc_list_up))
        intersect_genes_co_down=list(set(balf_list_down) & set(pbmc_list_down))
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
        print(ctype)
        print('n BALF: ' + str(n_balf))
        print('n PBMC: ' + str(n_pbmc))
        print('Union: '+str(U))
        print('Contra-Intersection: '+str(I_contra))
        print('Co-Intersection: '+str(I_co)+'\n')
        
        # Store codirectional intersecting genes to add to a dataframe
        int_list_up.append(intersect_genes_co_up)
        int_list_down.append(intersect_genes_co_down)
        
    # Create the intersection dataframe for upreg genes
    table_up = int_list_up
    int_df_up = pd.DataFrame(table_up)
    int_df_up = int_df_up.transpose()
    int_df_up.columns = overlap_celltypes
    
    # Print the dataframe to a .csv
    int_df_up.to_csv(res_dir+'STRINGENT_INT_DEG_UP_'+comparison+'_'+BALF_dataset+'_'+PBMC_dataset+'.csv',index=False)
    
    # Create the intersection dataframe for downreg genes
    table_down = int_list_down
    int_df_down = pd.DataFrame(table_down)
    int_df_down = int_df_down.transpose()
    int_df_down.columns = overlap_celltypes
    
    # Print the dataframe to a .csv
    int_df_down.to_csv(res_dir+'STRINGENT_INT_DEG_DOWN_'+comparison+'_'+BALF_dataset+'_'+PBMC_dataset+'.csv',index=False)
    
    df=df.T
    return(df)

    
def filter_genes(data):
    # Following DE remove genes we don't care about that could be artefactual
    # This will help to make more valid compartmental comparisons
    
    # Filter genes based on regex
    to_delete=data.filter(regex="^RP[SL]|^HB[^(P)]|^MALAT1$|^NEAT1$|^PPBP$|^MT-|^RNA18S|^RNA28S|^IGH|^IGL|^IGK|^TRAV|^TRBV|^TRDV|^TRGV", axis=0)
    
    data = data[~data.index.isin(to_delete.index)]

    return(data)

    
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

x = gets_degs(prefix = "Liao_", condition = 'Severe_Healthy Control_', celltypes=celltypes)   
y = matches2dataframe(prefix = "Liao_", condition = 'Severe_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Severe_vs_Healthy_Control_", sample_type = dataset, degs=x)

# Mild vs. Healthy Control
x = gets_degs(prefix = "Liao_", condition = 'Mild_Healthy Control_', celltypes=celltypes)   
y = matches2dataframe(prefix = "Liao_", condition = 'Mild_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Mild_vs_Healthy_Control_", sample_type = dataset, degs=x)

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

x = gets_degs(prefix = "Liao_",condition = 'Severe_Mild_', celltypes=celltypes)     
y = matches2dataframe(prefix = "Liao_", condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Severe_vs_Mild_", sample_type = dataset, degs=x)


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
            "Proliferating T"
]

x = gets_degs(prefix = "Wauters_", condition = 'Severe_Mild_', celltypes=celltypes)   
y = matches2dataframe(prefix = "Wauters_", condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Severe_vs_Mild_", sample_type = dataset, degs=x)

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

x = gets_degs(prefix = "Wilk_",condition = 'Severe_Healthy Control_', celltypes=celltypes)    
y = matches2dataframe(prefix = "Wilk_", condition = 'Severe_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Severe_vs_Healthy_Control_", sample_type = dataset, degs=x)


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

x = gets_degs(prefix = prefix,condition = 'Severe_Healthy Control_', celltypes=celltypes)    
y = matches2dataframe(prefix = prefix, condition = 'Severe_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Severe_vs_Healthy_Control_", sample_type = dataset, degs=x)

# Severe vs. Mild
x = gets_degs(prefix = prefix,condition = 'Severe_Mild_', celltypes=celltypes)               
y = matches2dataframe(prefix = prefix, condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Severe_vs_Mild_", sample_type = dataset, degs=x)

# Mild vs. Healthy Control
x = gets_degs(prefix = prefix,condition = 'Mild_Healthy Control_', celltypes=celltypes)               
y = matches2dataframe(prefix = prefix, condition = 'Mild_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Mild_vs_Healthy_Control_", sample_type = dataset, degs=x)


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

x = gets_degs(prefix = prefix,condition = 'Severe_Healthy Control_', celltypes=celltypes)             
y = matches2dataframe(prefix = prefix, condition = 'Severe_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Severe_vs_Healthy_Control_", sample_type = dataset, degs=x)

# Severe vs. Mild
x = gets_degs(prefix = prefix,condition = 'Severe_Mild_', celltypes=celltypes)        
y = matches2dataframe(prefix = prefix, condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Severe_vs_Mild_", sample_type = dataset, degs=x)

# Mild vs. Healthy Control
x = gets_degs(prefix = prefix,condition = 'Mild_Healthy Control_', celltypes=celltypes)        
y = matches2dataframe(prefix = prefix, condition = 'Mild_Healthy Control_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Mild_vs_Healthy_Control_", sample_type = dataset, degs=x)


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
x = gets_degs(prefix = prefix,condition = 'Severe_Mild_', celltypes=celltypes)        
y = matches2dataframe(prefix = prefix, condition = 'Severe_Mild_',celltypes=celltypes,ordered_celltypes = (celltypes),comparison = "Severe_vs_Mild_", sample_type = dataset, degs=x)

    
#-------------------------------------------------------------------------------------------------------------------------------#    
# All DEGs Compartmental Comparison Rose Plots [1]

########################## LIAO VS WILK ##########################
BALF_dataset="Liao_BALF"
PBMC_dataset="Wilk_PBMC"

# (1) Severe vs Healthy
comparison="Severe_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)


########################## LIAO VS LEE ##########################
BALF_dataset="Liao_BALF"
PBMC_dataset="Lee_PBMC"

# (2) Severe vs Healthy
comparison="Severe_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)

# (3) Severe vs. Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)

# (4) Mild vs. Healthy Control
comparison="Mild_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)


########################## LIAO VS ARUN ##########################
BALF_dataset="Liao_BALF"
PBMC_dataset="Arun_PBMC"

# (5) Severe vs Healthy
comparison="Severe_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)

# (6) Severe vs. Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)

# (7) Mild vs. Healthy Control
comparison="Mild_vs_Healthy_Control"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)


########################## LIAO VS SCHULTE ##########################
BALF_dataset="Liao_BALF"
PBMC_dataset="Schulte_PBMC"

# (8) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)

########################## WAUTERS VS LEE ##########################
BALF_dataset="Wauters_BALF"
PBMC_dataset="Lee_PBMC"

# (9) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)

########################## WAUTERS VS ARUN ##########################
BALF_dataset="Wauters_BALF"
PBMC_dataset="Arun_PBMC"

# (10) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)

########################## WAUTERS VS SCHULTE ##########################
BALF_dataset="Wauters_BALF"
PBMC_dataset="Schulte_PBMC"

# (11) Severe vs Mild
comparison="Severe_vs_Mild"
df=jaccard(BALF_dataset,PBMC_dataset,comparison)

