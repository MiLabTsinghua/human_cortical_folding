from itertools import chain

import anndata as ad
import itertools
import networkx as nx
import pandas as pd
import scanpy as sc
import scglue
import seaborn as sns
from matplotlib import rcParams
import os

scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)
os.chdir('/home/gaojie/workspace/Mida_collab/')


rna = sc.read_h5ad('Multi-omics_intergration/data/adata_RNA_final_input.h5ad')
atac = sc.read_h5ad('Multi-omics_intergration/data/adata_ATAC_final_input.h5ad')
guidance = nx.read_graphml("Multi-omics_intergration/data/guidance.graphml.gz")


scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep='X_pca_harmony',use_cell_type='celltype'
)
scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep='X_lsi_harmony',
)

guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()

glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": "glue"}
)

glue.save("Multi-omics_intergration/output/glue.dill")
rna.write('Multi-omics_intergration/output/adata_rna_output.h5ad')
atac.write('Multi-omics_intergration/output/adata_atac_output.h5ad')