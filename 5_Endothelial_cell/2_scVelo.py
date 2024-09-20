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
# os.chdir("/data/home/lixin/Project/SCC/Endothelial_cell/scVelo")
os.chdir("/5_Endothelial_cell/2_scVelo/")
# load sparse matrix:
X = io.mmread("Endothelial_cell_counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("Endothelial_cell_metadata.csv")

# load gene names:
with open("Endothelial_cell_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode'].values
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("Endothelial_cell_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['cell.type'], frameon=False, save=True)

# save dataset as anndata format
adata.write('Endothelial_cell.h5ad')

# reload dataset
adata = sc.read_h5ad('Endothelial_cell.h5ad')



scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2
adata = sc.read_h5ad('Endothelial_cell.h5ad')
# load loom files for spliced/unspliced matrices for each sample:

ldata = sc.read_h5ad("Endothelial_cell.loom.h5ad")
# ldata.write_loom("Endothelial_cell_loom.loom")
# ldata = scv.read("Endothelial_cell_loom.loom", cache=True)

ldata.obs.index = ldata.obs.barcode.tolist()
ldata.var_names_make_unique()

adata.obs.index = adata.obs.barcode.tolist()
adata.var_names_make_unique()
# merge matrices into the original adata object
# adata = scv.utils.merge(adata, ldata)
adata = scv.utils.merge(ldata, adata)

# plot umap to check
sc.pl.umap(adata, color='cell.type', frameon=False, legend_loc='on data', title='', save='Endothelial_cell_cluster.pdf')

#Part 2: Computing RNA velocity using scVelo
# scv.pl.proportions(adata, groupby='RNA_snn_res.0.3')

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

#Part 2.1: Visualize velocity fields
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='Endothelial_cell_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='cell.type', palette=["#644b8c", "#abbe75", "#bf9bba"], save='Endothelial_cell_embedding_grid.pdf', title='', scale=0.25)

plt.figure(figsize=(8, 8))
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['cell.type'], palette=["#644b8c", "#abbe75", "#bf9bba"], save='Endothelial_cell_embedding_stream.svg', title='', figsize=[3,4.2])
# "M2a"="#1F77B4FF","M2c"="#FF7F0EFF","CXCL1-3+ MDSCs"="#2CA02CFF","Monocyte"="#D62728FF","CD14+ DC"="#9467BDFF",
# "M2d"="#8C564BFF","CXCL9-11+ MDSCs"="#E377C2FF","CLEC9A+ DC"="#7F7F7FFF","CD1a+ CD1c+ DC"="#BCBD22FF","Cyclling M2"="#17BECFFF"
adata.obs["cell.type"]
adata.obs["cell.type"].unique()

#############visilization of gene expression along trajectory
#Similar to pseudotime, ‘latent time’ is computed from the dynamical model.
scv.tl.rank_velocity_genes(adata, groupby='cell.type', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

scv.tl.velocity_confidence(adata)

scv.tl.recover_dynamics(adata, n_jobs=32)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save="latent_time.svg")

#输出伪时间
adata.obs['latent_time']
adata.obs['velocity_pseudotime']

latent_time = adata.obs['latent_time'].to_frame()
print(type(latent_time))
velocity_pseudotime = adata.obs['velocity_pseudotime'].to_frame()
dt = pd.concat([latent_time, velocity_pseudotime], axis=1)
dt.to_csv("Endothelial_latent_time.txt", sep='\t')
