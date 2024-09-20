import loompy as lp
# print(loompy.__version__)
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import scvelo as scv
import cellrank as cr
import IPython
import matplotlib
import matplotlib.pyplot as plt

# set work path
os.chdir("/4_Myeloid/3_scVelo/")
# load sparse matrix:
X = io.mmread("Myeloid_cell_counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("Myeloid_cell_metadata.csv")

# load gene names:
with open("Myeloid_cell_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode'].values
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("Myeloid_cell_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['cell.type'], frameon=False, save=True)

# save dataset as anndata format
adata.write('Myeloid_cell.h5ad')

# reload dataset
adata = sc.read_h5ad('Myeloid_cell.h5ad')



scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('Myeloid_cell.h5ad')
# load loom files for spliced/unspliced matrices for each sample:

ldata = sc.read_h5ad("Myeloid_cell.loom.h5ad")
# ldata.write_loom("Myeloid_cell_loom.loom")
# ldata = scv.read("Myeloid_cell_loom.loom", cache=True)

ldata.obs.index = ldata.obs.barcode.tolist()
ldata.var_names_make_unique()

adata.obs.index = adata.obs.barcode.tolist()
adata.var_names_make_unique()
# merge matrices into the original adata object
# adata = scv.utils.merge(adata, ldata)
adata = scv.utils.merge(ldata, adata)

# plot umap to check
sc.pl.umap(adata, color='cell.type', frameon=False, legend_loc='on data', title='', save='Myeloid_cell_cluster.pdf')

#Part 2: Computing RNA velocity using scVelo
# scv.pl.proportions(adata, groupby='RNA_snn_res.0.3')

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

#Part 2.1: Visualize velocity fields
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='Myeloid_cell_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='cell.type', palette=["#BA94C7","#FCB8B9","#ECBA7D","#EEDFCA","#F6F195","#A7BBD6"], save='Myeloid_cell_embedding_grid.pdf', title='', scale=0.25)

plt.figure(figsize=(8, 8))
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['cell.type'], palette=["#BA94C7","#FCB8B9","#ECBA7D","#EEDFCA","#F6F195","#A7BBD6"], save='Myeloid_cell_embedding_stream.svg', title='')
# "M2a"="#1F77B4FF","M2c"="#FF7F0EFF","CXCL1-3+ MDSCs"="#2CA02CFF","Monocyte"="#D62728FF","CD14+ DC"="#9467BDFF",
# "M2d"="#8C564BFF","CXCL9-11+ MDSCs"="#E377C2FF","CLEC9A+ DC"="#7F7F7FFF","CD1a+ CD1c+ DC"="#BCBD22FF","Cyclling M2"="#17BECFFF"
adata.obs["cell.type"]
adata.obs["cell.type"].unique()
