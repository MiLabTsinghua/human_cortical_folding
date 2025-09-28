import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import anndata as ad
import celloracle as co
import scipy
from matplotlib import cm
from matplotlib.colors import Normalize, to_hex
import celloracle_fun
co.__version__

plt.rcParams['figure.figsize'] = [4, 4]
plt.rcParams["savefig.dpi"] = 300
os.chdir('/home/gaojie/workspace/Mida_collab/TF_perturbation')


oracle = co.load_hdf5("output/GW11.celltype_sample.celloracle.oracle")
pseudoage = pd.read_csv('../sample_data/RNA_raw/pseudoage_folding.vRG.tsv',sep='\t')
oracle.adata = oracle.adata[oracle.adata.obs['celltype']=='vRG',:] ###subset to vRG
oracle.adata.obs['pseudoage'] = pseudoage.loc[oracle.adata.obs_names,'PC1'] ###adding pseudoage info


links = co.load_hdf5("output/links.GW11.celltype_sample.celloracle.links")
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
links.cluster = ['vRG_Adjacent',
 'vRG_Distant',
 'vRG_Sulcus']

# Calculate network scores.     
links.get_network_score()   
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=20,
                              use_cluster_specific_TFdict=True)

oracle_raw = oracle

all_TFs_df = pd.read_csv('../CellOracle/output/GSEA_res_df_exp0.04pct0.05size20.csv',index_col=0)
all_TFs_df = all_TFs_df[all_TFs_df['padj.2']<=0.05]
all_TFs = all_TFs_df['tf'].unique().tolist()
print(all_TFs,flush=True)
print(len(all_TFs),flush=True)


gradient = celloracle_fun.prepare_gradient(oracle = oracle_raw,pseudotime_key = 'pseudoage') ###see celloracle_fun.py

###test all TF effect
for TF in all_TFs:
    print(TF)
    try:
        oracle = celloracle_fun.do_perturb(oracle=oracle_raw,TF=TF)
        draw_TF_perturabation(oracle=oracle,TF=TF,min_mass=18) ###draw perturbation UMAP 
        dev = celloracle_fun.calculate_PS(oracle = oracle,gradient = gradient,TF = TF) ###save PS score csv
    except Exception as e:
        print(Exception)
        continue

