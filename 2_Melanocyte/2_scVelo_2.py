import loompy as lp
print(lp.__version__)
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

# set work path
os.chdir("/2_Melanocyte/2_scVelo")
# load sparse matrix:
X = io.mmread("Melanocyte_counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("Melanocyte_metadata.csv")

# load gene names:
with open("Melanocyte_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode'].values
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("Melanocyte_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['RNA_snn_res.0.3'], frameon=False, save=True)

# save dataset as anndata format
adata.write('Melanocyte.h5ad')

# reload dataset
adata = sc.read_h5ad('Melanocyte.h5ad')



scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('Melanocyte.h5ad')
# load loom files for spliced/unspliced matrices for each sample:

ldata = sc.read_h5ad("Melanocyte_loom.h5ad")
# ldata.write_loom("Melanocyte_loom.loom")
# ldata = scv.read("Melanocyte_loom.loom", cache=True)

ldata.obs.index = ldata.obs.barcode.tolist()
ldata.var_names_make_unique()

adata.obs.index = adata.obs.barcode.tolist()
adata.var_names_make_unique()
# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

# plot umap to check
sc.pl.umap(adata, color='RNA_snn_res.0.3', frameon=False, legend_loc='on data', title='', save='Melanocyte_cluster.pdf')

#Part 2: Computing RNA velocity using scVelo
# scv.pl.proportions(adata, groupby='RNA_snn_res.0.3')

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

#Part 2.1: Visualize velocity fields
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='Melanocyte_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='RNA_snn_res.0.3', palette=["#db6657","#977dac","#4876b4","#b1ccea"], save='Melanocyte_embedding_grid.pdf', title='', scale=0.25)

scv.pl.velocity_embedding_stream(adata, basis='umap', color=['RNA_snn_res.0.3'], palette=["#db6657","#977dac","#b1ccea","#4876b4"], save='Melanocyte_embedding_stream.svg', title='', figsize=[6,4.7])
adata.obs["RNA_snn_res.0.3"].unique()



