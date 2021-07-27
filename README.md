# 2021_iScience_Overholt_Krog_COVID

# About

This repository contains code related to the manuscript titled  "Dissecting the common and compartment-specific features of COVID-19 severity in the lung and periphery with single-cell resolution" by K.J. Overholt, J.R. Krog, I. Zanoni, and B.D. Bryson in *iScience* (2021).

DOI: https://doi.org/10.1016/j.isci.2021.102738

## Contents
1. **R** - Contains R scripts
2. **Python** - Contains Python scripts

## scRNA-seq Datasets
See manuscript for GEO (https://www.ncbi.nlm.nih.gov/geo) and EGA (https://ega-archive.org) download details.



## Order for running scripts
1.	main_int_harmony.R *
2.	main_int_clustering.R
3.	subcluster_myeloid.R
4.	subcluster_NKT.R
5.	subcluster_PCBPro.R
6.	subcluster_dendritic.R
7.	recombine.R
8.	de_intrasev.R
9.	de_\<cohort>\_rna.R (run for all cohorts)
10.	de_\<cohort>\_rna_donor.R (run for all cohorts)
11.	stringent_degs.py
12.	all_degs_plots.py
13.	surface_degs_plots.py \**
14.	intrasev_plots.py
15.	gsea_\<cohort>\.R (run for all cohorts) †
16.	gsea_intrasev.R †
17.	gsea_plots.py
18.	gsea_intrasev_plots.py
19.	enrichr_intrasev.R
20.	enrichr_heatmap.py
21.	nichenet_loop.R
22.	printr.R (run for each cohort desired)
23.	ligands_broad.py
24.	ligands_targeted.py
  
 ## Necessary accessory files
  
 \* meta.txt at https://github.com/zhangzlab/covid_balf/blob/master/meta.txt
 
 \** human_surfaceproteins.xlsx from "CSPA validated surfaceome proteins" at https://wlab.ethz.ch/cspa/#downloads
  
 †  h.all.v7.1.symbols.gmt from "H: hallmark gene sets" at https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
