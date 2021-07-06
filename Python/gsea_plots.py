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

IFN_dir = wd + '/GSEA/IFN_Heatmaps/'


# ### Plotting function for Dotplot

def dotplot(pw_by_cells,p_vals,size,option,title):
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
    fig.savefig(res_dir+'gsea_dotplot_'+title+'_'+option+'.pdf')


# #### Plotting function for Heatmap

def heatmap(pw_by_cells,title,option):
    fig= figure(num=None, figsize=(9, 8), dpi=300)
    ax  = fig.add_subplot(111)

    pw_by_cells = pw_by_cells[pw_by_cells.columns].astype(float) # Change the values from 'object' to 'float'
    pw_by_cells = pw_by_cells.reindex(index=pw_by_cells.index[::-1]) # Reverse the order of the rows so that it matches the dotplot
    # Also try cmap RdBu_r

    sns.heatmap(pw_by_cells,square=False,xticklabels=True,yticklabels=True,cmap="coolwarm",cbar_kws={'label': 'Normalized Enrichment Score'},fmt='',center=0)
    fig.tight_layout()
    fig.savefig(res_dir+'gsea_heatmap_'+title+'_'+option+'.pdf')


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


def interferons(filenames,celltypes,option,title):

    J_df=pd.DataFrame(columns=['Jaccard'],index=celltypes)
    unique_alpha_masterlist=[]
    unique_gamma_masterlist=[]

    for filename, cell in zip(filenames, celltypes):

        df=pd.read_csv(filename,sep='\t')  
        
        # Find leading edge of the IFN_ALPHA and IFN_GAMMA pathways
        # Store the Jaccard index and unique elements of both pathways
        J,unique_alpha,unique_gamma=leading_edge(df)
        J_df.loc[cell]['Jaccard']=J
        
        # Extend (not append) unique_alpha and unique_gamma to their own lists
        unique_alpha_masterlist.extend(unique_alpha)
        unique_gamma_masterlist.extend(unique_gamma)
   
    # Make the lists unique
    unique_alpha_masterlist=list(set(unique_alpha_masterlist))
    unique_gamma_masterlist=list(set(unique_gamma_masterlist))
    
    # Generate unique_alpha df and unique_gamma df
    alpha_df=pd.DataFrame(index=celltypes,columns=unique_alpha_masterlist)
    gamma_df=pd.DataFrame(index=celltypes,columns=unique_gamma_masterlist)
    
    for filename,cell in zip(filenames,celltypes):
        df=pd.read_csv(filename,sep='\t') 
        J,unique_alpha,unique_gamma=leading_edge(df)
        for gene_alpha in unique_alpha:
            alpha_df.loc[cell][gene_alpha]=1
        for gene_gamma in unique_gamma:
            gamma_df.loc[cell][gene_gamma]=1
            
    ledge_heatmap(alpha_df,gamma_df,J_df,option,title)

    return()


def leading_edge(df):
    
    df=df.set_index('pathway')
    
    if ('HALLMARK_INTERFERON_ALPHA_RESPONSE' in df.index and 'HALLMARK_INTERFERON_GAMMA_RESPONSE' in df.index):
        df=df.loc[['HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE']]

        if (df.loc['HALLMARK_INTERFERON_ALPHA_RESPONSE']['padj']<0.05 and df.loc['HALLMARK_INTERFERON_GAMMA_RESPONSE']['padj']<0.05):
            ledge_alpha=df.loc['HALLMARK_INTERFERON_ALPHA_RESPONSE']['leadingEdge'].split() # Split is needed to convert string to list
            ledge_gamma=df.loc['HALLMARK_INTERFERON_GAMMA_RESPONSE']['leadingEdge'].split()

            intersect=list(set(ledge_alpha) & set(ledge_gamma))

            I=len(intersect)
            U=len(ledge_alpha)+len(ledge_gamma)-I
            J=I/U

            unique_alpha=list(set(ledge_alpha) - set(intersect))
            unique_gamma=list(set(ledge_gamma) - set(intersect))

        else:
            J=0
            unique_alpha=[]
            unique_gamma=[]
    else:
            J=0
            unique_alpha=[]
            unique_gamma=[]
            
    return(J,unique_alpha,unique_gamma)

def ledge_heatmap(alpha_df,gamma_df,J_df,option,title):

    # Fill NAs with zeros
    alpha_df=alpha_df.fillna(0).T
    gamma_df=gamma_df.fillna(0).T
    J_df=J_df.fillna(0).T
    
    # Set scale for the heatmaps
    scale=0.5

    fig= figure(figsize=(scale*alpha_df.shape[1],scale*alpha_df.shape[0]))
    ax  = fig.add_subplot(111)
    sns.heatmap(alpha_df,square=True,xticklabels=True,yticklabels=True,cmap="YlGnBu",cbar_kws={'label': 'Binary'},fmt='',vmin=0,vmax=1)
    fig.tight_layout()
    fig.savefig(IFN_dir+'gsea_IFNA_'+title+'_'+option+'.pdf')
    
    fig= figure(figsize=(scale*gamma_df.shape[1],scale*gamma_df.shape[0]))
    ax  = fig.add_subplot(111)
    sns.heatmap(gamma_df,square=True,xticklabels=True,yticklabels=True,cmap="PuRd",cbar_kws={'label': 'Binary'},fmt='',vmin=0,vmax=1)
    fig.tight_layout()
    fig.savefig(IFN_dir+'gsea_IFNG_'+title+'_'+option+'.pdf')
    
    fig= figure(figsize=(scale*J_df.shape[1],scale*10))
    ax  = fig.add_subplot(111)
    sns.heatmap(J_df,square=True,xticklabels=True,yticklabels=True,cmap="BuPu",cbar_kws={'label': 'Jaccard Index',"shrink": 1.0},fmt='',vmin=0,vmax=1)
    fig.tight_layout()
    fig.savefig(IFN_dir+'gsea_IFN_JAC_'+title+'_'+option+'.pdf')
    
    J_df.to_csv(IFN_dir+'IFN_JAC_'+title+'_'+option+'.csv')
    

def jaccard(p_vals_BALF,p_vals_PBMC,pw_by_cells_BALF,pw_by_cells_PBMC):
    
    pw_by_cells_BALF[p_vals_BALF.values<(-log10(0.05))] = 0
    pw_by_cells_PBMC[p_vals_PBMC.values<(-log10(0.05))] = 0
    
    cells_by_pw_BALF = pw_by_cells_BALF.T
    cells_by_pw_PBMC = pw_by_cells_PBMC.T
    
    # Find overlapping cell types in BALF and PBMC
    overlap_celltypes=list(set(cells_by_pw_BALF.index) & set(cells_by_pw_PBMC.index))

    df = pd.DataFrame(columns = overlap_celltypes,index=['Frac BALF Only','Frac PBMC Only','Jaccard Contradirectional','Jaccard Codirectional'])

    for ctype in overlap_celltypes:
        # Filter for only nonzero columns in each BALF cell type, count how many
        # Create separate lists for upreg and downreg pathways
        foo=cells_by_pw_BALF.loc[ctype]
        balf_list_up=foo[foo>0].index
        balf_list_down=foo[foo<0].index
        n_balf_up=len(balf_list_up)
        n_balf_down=len(balf_list_down)
        n_balf=n_balf_up+n_balf_down

        # Filter for only nonzero columns in each PBMC cell type, count how many
        # Create separate lists for upreg and downreg pathways
        foo=cells_by_pw_PBMC.loc[ctype]
        pbmc_list_up=foo[foo>0].index
        pbmc_list_down=foo[foo<0].index
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
        print(ctype)
        print('n BALF: ' + str(n_balf))
        print('n PBMC: ' + str(n_pbmc))
        print('Union: '+str(U))
        print('Contra-Intersection: '+str(I_contra))
        print('Co-Intersection: '+str(I_co)+'\n')

    df=df.T
    return(df)


# ## Nightingale Rose Plot
# See https://medium.com/@abhishekdas.nitc/nightingale-plots-in-python-using-plotly-da42bc18d15d
# https://plotly.com/python/polar-chart/#polar-bar-chart
# https://plotly.com/python/static-image-export/

def roseplot(df,sets,comparison):

    fig = go.Figure()

    color_J_co='rgb(253,231,37)'
    color_J_contra='rgb(109,205,89)'
    color_P='rgb(31,158,137)'
    color_B='rgb(68,1,84)'

#     color_J_co='rgb(253,231,37)'
#     color_J_contra='rgb(214, 255, 225)'
#     color_P='rgb(214, 238, 255)'
#     color_B='rgb(239, 214, 255)'

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
        title=sets+' '+comparison+' GSEA Pathways',
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
    plotly.offline.plot(fig, filename=res_dir+'rose_GSEA_'+sets+'_'+comparison+'.html')
    fig.write_image(res_dir+'rose_GSEA_'+sets+'_'+comparison+'.pdf') # Must have Kaleido installed


# Function to prepare for dotplot, heatmap, and interferon plots
# Calls dotplot, interferons, and reader
def prep_plots(option,celltypes,dataset,size):
    filenames=[]
    # Dataset should be of the form "Lee_"
    # size is (x,y)

    # Severe vs. Mild option
    if option=='SM':
        for cell in celltypes:
            filenames.append(data_dir+dataset+cell+" Severe Mild fgsea.tsv")

    # Mild vs. Control option      
    if option=='MC':
        for cell in celltypes:
            filenames.append(data_dir+dataset+cell+" Mild Healthy Control fgsea.tsv")

    # Severe vs. Control option      
    if option=='SC':
        for cell in celltypes:
            filenames.append(data_dir+dataset+cell+" Severe Healthy Control fgsea.tsv")        
    cols_rename=celltypes

    pw_by_cells,p_vals=reader(filenames,celltypes,cols_rename)
    
    # Make the desired plots
    #heatmap(pw_by_cells_BALF,title,option)
    dotplot(pw_by_cells,p_vals,size,option,dataset)
    interferons(filenames,celltypes,option,dataset)
    
    return(p_vals,pw_by_cells)


def runner(option,celltypes_BALF,celltypes_PBMC,BALF_dataset,PBMC_dataset):
    
    p_vals_BALF,pw_by_cells_BALF=prep_plots(option,celltypes_BALF,BALF_dataset,(7,10))
    p_vals_PBMC,pw_by_cells_PBMC=prep_plots(option,celltypes_PBMC,PBMC_dataset,(7,10))
    df=jaccard(p_vals_BALF,p_vals_PBMC,pw_by_cells_BALF,pw_by_cells_PBMC)
    
    title=BALF_dataset+PBMC_dataset
    roseplot(df,title,option)

#-------------------------------------------------------------------------------------------------------------------------------#    

# Define celltypes
celltypes_liao_SC = [
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

celltypes_liao_MC=celltypes_liao_SC

celltypes_liao_SM = [
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

celltypes_wilk = [
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

celltypes_arun = [
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

celltypes_lee = [
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

celltypes_schulte = [
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

celltypes_wauters = [
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

#-------------------------------------------------------------------------------------------------------------------------------#  
# Cross-Compartmental comparison rose plots, dotplots, heatmaps, interferon plots [1]

########################## LIAO VS WILK ##########################
BALF_dataset="Liao_"
PBMC_dataset="Wilk_"

# (1) Severe vs Healthy
option='SC'
runner(option,celltypes_liao_SC,celltypes_wilk,BALF_dataset,PBMC_dataset)

########################## LIAO VS LEE ##########################
BALF_dataset="Liao_"
PBMC_dataset="Lee_"

# (2) Severe vs Healthy
option='SC'
runner(option,celltypes_liao_SC,celltypes_lee,BALF_dataset,PBMC_dataset)

# (3) Severe vs. Mild
option='SM'
runner(option,celltypes_liao_SM,celltypes_lee,BALF_dataset,PBMC_dataset)

# (4) Mild vs. Healthy Control
option='MC'
runner(option,celltypes_liao_MC,celltypes_lee,BALF_dataset,PBMC_dataset)


########################## LIAO VS ARUN ##########################
BALF_dataset="Liao_"
PBMC_dataset="Arun_"

# (5) Severe vs Healthy
option='SC'
runner(option,celltypes_liao_SC,celltypes_arun,BALF_dataset,PBMC_dataset)

# (6) Severe vs. Mild
option='SM'
runner(option,celltypes_liao_SM,celltypes_arun,BALF_dataset,PBMC_dataset)

# (7) Mild vs. Healthy Control
option='MC'
runner(option,celltypes_liao_MC,celltypes_arun,BALF_dataset,PBMC_dataset)


########################## LIAO VS SCHULTE ##########################
BALF_dataset="Liao_"
PBMC_dataset="Schulte_"

# (8) Severe vs Mild
option='SM'
runner(option,celltypes_liao_SM,celltypes_schulte,BALF_dataset,PBMC_dataset)

########################## WAUTERS VS LEE ##########################
BALF_dataset="Wauters_"
PBMC_dataset="Lee_"

# (9) Severe vs Mild
option='SM'
runner(option,celltypes_wauters,celltypes_lee,BALF_dataset,PBMC_dataset)

########################## WAUTERS VS ARUN ##########################
BALF_dataset="Wauters_"
PBMC_dataset="Arun_"

# (10) Severe vs Mild
option='SM'
runner(option,celltypes_wauters,celltypes_arun,BALF_dataset,PBMC_dataset)

########################## WAUTERS VS SCHULTE ##########################
BALF_dataset="Wauters_"
PBMC_dataset="Schulte_"

# (11) Severe vs Mild
option='SM'
runner(option,celltypes_wauters,celltypes_schulte,BALF_dataset,PBMC_dataset)


#-------------------------------------------------------------------------------------------------------------------------------# 
# Cross-Compartmental comparison rose plots, dotplots, heatmaps, interferon plots [6]


########################## LIAO VS WAUTERS ##########################
dataset1="Liao_"
dataset2="Wauters_"

# (1) Severe vs Mild
option="SM"
runner(option,celltypes_liao_SM,celltypes_wauters,dataset1,dataset2)

########################## WILK VS LEE ##########################
dataset1="Wilk_"
dataset2="Lee_"

# (2) Severe vs Healthy
option="SC"
runner(option,celltypes_wilk,celltypes_lee,dataset1,dataset2)

########################## WILK VS ARUN ##########################
dataset1="Wilk_"
dataset2="Arun_"

# (3) Severe vs Healthy
option="SC"
runner(option,celltypes_wilk,celltypes_arun,dataset1,dataset2)

########################## LEE VS ARUN ##########################
dataset1="Lee_"
dataset2="Arun_"

# (4) Severe vs Healthy
option="SC"
runner(option,celltypes_lee,celltypes_arun,dataset1,dataset2)

# (5) Severe vs Mild
option="SM"
runner(option,celltypes_lee,celltypes_arun,dataset1,dataset2)

# (6) Mild vs Healthy
option="MC"
runner(option,celltypes_lee,celltypes_arun,dataset1,dataset2)

########################## LEE VS SCHULTE ##########################
dataset1="Lee_"
dataset2="Schulte_"

# (7) Severe vs Mild
option="SM"
runner(option,celltypes_lee,celltypes_schulte,dataset1,dataset2)

########################## ARUN VS SCHULTE ##########################
dataset1="Arun_"
dataset2="Schulte_"

# (8) Severe vs Mild
option="SM"
runner(option,celltypes_arun,celltypes_schulte,dataset1,dataset2)

