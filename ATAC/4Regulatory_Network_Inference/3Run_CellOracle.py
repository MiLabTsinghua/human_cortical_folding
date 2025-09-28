import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import anndata as ad
import celloracle as co
co.__version__

os.chdir('/home/gaojie/workspace/Mida_collab/')
adata = sc.read_h5ad("CellOracle/input_/scRNA_input_celltype.h5ad")
print(adata.var['highly_variable'].sum())

adata.raw = adata 
adata = adata[:, adata.var.highly_variable]

print(f"Cell number is :{adata.shape[0]}")
print(f"Gene number is :{adata.shape[1]}")

print("Metadata columns :", list(adata.obs.columns))
print("Dimensional reduction: ", list(adata.obsm.keys()))

sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

adata.X = adata.layers["counts"].copy()

base_GRN = pd.read_parquet("CellOracle/input/base_GRN_dataframe.parquet")

### cluster = celltype
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="celltype",
                                   embedding_name="X_umap")
# You can load TF info dataframe with the following code.
oracle.import_TF_data(TF_info_matrix=base_GRN)
oracle.perform_PCA()

# Select important PCs
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
print(n_comps)
n_comps = min(n_comps, 50)
n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=20)
oracle.to_hdf5("CellOracle/output/GW11.celltype.celloracle.oracle") ### must end with .celloracle.oracle

links = oracle.get_links(cluster_name_for_GRN_unit="celltype", alpha=10,
                         verbose_level=10,n_jobs=20)

links.to_hdf5(file_path="CellOracle/output/links.GW11.celltype.celloracle.links")### must end with .celloracle.links