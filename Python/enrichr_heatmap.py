import numpy as np
import pandas as pd
import scipy as sp
import statistics
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import preprocessing
import os

# Get the directories for reading/saving files
wd = os.getcwd()

enrich_dir = wd + '/GO/EnrichR_Results/'
res_dir = wd + '/GO/EnrichR_Plots/'

if not os.path.exists(res_dir):
    os.mkdir(res_dir)


# # Run EnrichR on OVERLAPPING DEGs between compartments

prefix_list = [
    "Severe_Mild",
    "Severe_Healthy Control",
    "Mild_Healthy Control"]

file_types = ["GO_Mol_","GO_BIO_","TF_"]

directions = ["UP_","DOWN_"]

for ftype in file_types:
    for direction in directions:
        for prefix in prefix_list:
            
            file_head=ftype+direction+prefix
            
            # Grab all files starting with a given file prefix
            fnames = [x for x in os.listdir(enrich_dir) if x.startswith(file_head)] 

            paths = set()
            celltypes=list()
            
            for fname in fnames:
                # Get the name of the celltype matching the file
                celltype=fname.replace(file_head, '')
                celltype=celltype.replace('.txt', '')
                celltypes.append(celltype)
                
                # Get the important GO terms
                data = pd.read_csv(enrich_dir+fname, sep="\t",index_col=1)
                top_data = data.loc[(data['Adjusted.P.value']<.05)].index 
                top_data = set(top_data)
                paths.update(top_data)
                
            paths = list(paths)
            
            # Create empty dataframe
            df = pd.DataFrame(index = celltypes,columns = paths)
            
            # Loop through the files and fill the dataframe
            for celltype,fname in zip(celltypes,fnames):
                for path in paths:
                    data = pd.read_csv(enrich_dir+fname, sep="\t",index_col=1)
                    if path in data.index and (data.loc[path]['Adjusted.P.value']<.05):
                        df.loc[celltype,path] = data.loc[path,"Combined.Score"]
                        
            # If you want to normalize 
#             normalized_df=(df.T-df.T.min())/((df.T.max()+1)-df.T.min())
#             normalized_df = normalized_df.T.fillna(0)
            
            # If not normalizing
            normalized_df=df.fillna(0)
        
            # Either step
            normalized_df = normalized_df.rename(columns=lambda s: s.rsplit(' ',1)[0])
            
            if (not normalized_df.empty):
            
                # sns.set(font_scale=1) 
                if direction == "UP_":
                    fig = plt.figure(figsize=(150,15))
                    g= sns.heatmap(normalized_df, annot=False,  cmap= 'viridis', cbar_kws={'label': "Normalized Enrichr Combined Score"},vmin=0, vmax=2000,square=True)
#                     g.figure.axes[-1].yaxis.label.set_size(12)
#                     g.set_yticklabels(g.get_ymajorticklabels(), fontsize = 12)
                    plt.title(ftype+direction+prefix,fontsize= 12,weight='bold')
                    plt.yticks(rotation=0)
                    plt.tight_layout()
                    plt.savefig(res_dir+ftype+direction+prefix+'.pdf')
                    plt.close("all")

                if direction == "DOWN_":
                    fig = plt.figure(figsize=(150,15))
                    g= sns.heatmap(normalized_df, annot=False,  cmap= 'viridis', cbar_kws={'label': "Normalized Enrichr Combined Score"},vmin=0, vmax=2000,square=True)
#                     g.figure.axes[-1].yaxis.label.set_size(12)
#                     g.set_yticklabels(g.get_ymajorticklabels(), fontsize = 12)
                    plt.title(ftype+direction+prefix,fontsize= 12,weight='bold')
                    plt.yticks(rotation=0)
                    plt.tight_layout()
                    plt.savefig(res_dir+ftype+direction+prefix+'.pdf')
                    plt.close("all")
            else:
                    print('This dataframe was EMPTY')
                
            print("Done with",ftype+direction+prefix)
            
print("Finished!")

