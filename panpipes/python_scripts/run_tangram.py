'''
Run Tangram
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg 
import muon as mu

import os
import argparse
import sys
import logging
import json


L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("testing logger works")



sc.settings.verbosity = 3

parser = argparse.ArgumentParser()

parser.add_argument("--input_spatial",
                    help="path to mudata of spatial transcriptomics data")
parser.add_argument("--input_singlecell",
                    help="path to mudata of single-cell reference data")
parser.add_argument("--figdir",
                    default="./figures/Tangram",
                    help="path to save the figures to")
parser.add_argument("--output_dir",
                    default="./tangram.output",
                    help="path to save the output to")

# parameters for feature selection: 
parser.add_argument("--gene_list",
                    default=None,
                    help="path to a csv containing a list of genes to use for the deconvolution. csv-file needs to contain a header")
parser.add_argument("--labels_key_rank_genes",
                    default=None,
                    help="which column in .obs of the reference to use for the 'groupby' parameter of sc.tl.rank_genes_groups()")
parser.add_argument("--n_genes_rank",
                    default=100,
                    help="how many top genes to select of each 'groupby' group")
parser.add_argument("--layer_rank_genes",
                    default=None,
                    help="which layer to use of the reference for sc.tl.rank_genes_groups(). if None, uses .X")
parser.add_argument("--method_rank_genes",
                    default=None,
                    help="which test method to use. one of: ['logreg', 't-test', 'wilcoxon', 't-test_overestim_var']")
parser.add_argument("--corr_method_rank_genes",
                    default="benjamini-hochberg",
                    help="which p-value correction method to use. one of: ['benjamini-hochberg', 'bonferroni']. Used only for 't-test', 't-test_overestim_var', and 'wilcoxon'")

# model parameters
parser.add_argument("--labels_key_model",
                    default=None,
                    help="cell type key in the reference .obs")
parser.add_argument("--num_epochs",
                    default=1000,
                    help="Number of epochs for tangram.mapping_utils.map_cells_to_space()")
parser.add_argument("--device",
                    default="cpu",
                    help="Device to use for deconvolution")
parser.add_argument("--kwargs",
                    default='{}',
                    type=str,
                    help="Parameters for tangram.mapping_utils.map_cells_to_space()")



args, opt = parser.parse_known_args()

L.info("running with args:")
L.debug(args)

if isinstance(args.kwargs, str): 
	args.kwargs = json.loads(args.kwargs)

figdir = args.figdir
if not os.path.exists(figdir):
    os.mkdir(figdir)
sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

output_dir = args.output_dir
if not os.path.exists(output_dir):
    os.mkdir(output_dir)



#1. read in the data
#spatial: 
mdata_spatial = mu.read(args.input_spatial)
adata_st = mdata_spatial.mod['spatial']
#single-cell: 
mdata_singlecell = mu.read(args.input_singlecell)
adata_sc = mdata_singlecell.mod['rna']


#2. Perform gene selection:
if args.gene_list is not None: # read in csv and create list 
    markers = pd.read_csv(args.gene_list, header = 0)
    markers = list(markers.iloc[:, 0])

else: # perform feature selection using sc.tl.rank_genes_groups()
    sc.tl.rank_genes_groups(adata_sc, groupby=args.labels_key_rank_genes, layer=args.layer_rank_genes, method=args.method_rank_genes,corr_method = args.corr_method_rank_genes)
    sc.pl.rank_genes_groups(adata_sc, show = False, save = ".png")
    markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:int(args.n_genes_rank), :]
    markers_df.to_csv(output_dir + "/rank_genes_groups.csv")
    markers = list(np.unique(markers_df.melt().value.values))

# "Preprocess" anndatas
tg.pp_adatas(adata_sc=adata_sc, adata_sp=adata_st, genes=markers)

# 3. Run tangram
adata_results = tg.mapping_utils.map_cells_to_space(
        adata_sc=adata_sc, adata_sp=adata_st, num_epochs=int(args.num_epochs), device=args.device, **args.kwargs
    )

# 3. Extract and plot results 
tg.project_cell_annotations(adata_results, adata_st, annotation=args.labels_key_model)

annotation_list = list(pd.unique(adata_sc.obs[args.labels_key_model]))
df = adata_st.obsm["tangram_ct_pred"][annotation_list]
tg.construct_obs_plot(df, adata_st, perc=0.05)
if "spatial" in adata_st.uns: 
	sc.pl.spatial(adata_st, color=annotation_list, cmap="viridis", show=False, frameon=False, ncols=3, save = "_tangram_ct_pred.png")
else: 
	sc.pl.spatial(adata_st, color=annotation_list, cmap="viridis", show=False, frameon=False, ncols=3, save = "_tangram_ct_pred.png",spot_size=0.5)


mdata_singlecell_results = mu.MuData({"rna": adata_sc})
mdata_spatial_results = mu.MuData({"spatial": adata_st})

mdata_singlecell_results.write(output_dir+"/Tangram_screference_output.h5mu")
mdata_spatial_results.write(output_dir+"/Tangram_spatial_output.h5mu")

L.info("Done")
