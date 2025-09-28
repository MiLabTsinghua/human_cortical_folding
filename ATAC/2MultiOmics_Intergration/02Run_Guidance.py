import scanpy as sc
import pandas as pd
import os
from scipy.io import mmread
from scipy.sparse import coo_matrix
import scglue
import numpy as np
from matplotlib import rcParams
import networkx as nx


os.chdir('/home/gaojie/workspace/Mida_collab/')

adata_atac = sc.read_h5ad('Multi-omics_intergration/data/adata_ATAC_final_input.h5ad')

adata_rna = sc.read_h5ad('Multi-omics_intergration/data/adata_RNA_final_input.h5ad')

guidance = scglue.genomics.rna_anchored_guidance_graph(adata_rna,adata_atac)
scglue.graph.check_graph(guidance, [adata_rna, adata_atac])

nx.write_graphml(guidance, "Multi-omics_intergration/data/guidance.graphml.gz")
