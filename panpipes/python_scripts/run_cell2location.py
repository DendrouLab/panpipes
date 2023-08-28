'''
Run Cell2Location
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import cell2location as c2l
import scanpy as sc
import pandas as pd
import muon as mu

import os
import argparse
import sys
import logging

from panpipes.funcs.plotting import cell2loc_plot_QC_reference
from panpipes.funcs.plotting import cell2loc_plot_QC_reconstr
from panpipes.funcs.plotting import cell2loc_plot_history
from panpipes.funcs.scmethods import cell2loc_filter_genes


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
                    default="./figures/Cell2Location",
                    help="path to save the figures to")
parser.add_argument("--output_dir",
                    default="./cell2location.output",
                    help="path to save the figures to")
parser.add_argument("--save_models",
                    default=False,
                    help="whether to save the reference & spatial mapping models")


# parameters for feature selection: 
parser.add_argument("--gene_list",
                    default="None",
                    help="path to a csv containing a list of genes to use for the deconvolution")
parser.add_argument("--remove_mt",
                    default=True,
                    help="whether to remove MT genes before deconvolution")
parser.add_argument("--cell_count_cutoff",
                    default=15,
                    help="all genes detected in less than cell_count_cutoff cells will be excluded")
parser.add_argument("--cell_percentage_cutoff2",
                    default=0.05,
                    help="all genes detected in at least this percentage of cells will be included")
parser.add_argument("--nonz_mean_cutoff",
                    default=1.12,
                    help="genes detected in the number of cells between the above mentioned cutoffs are selected only when their average expression in non-zero cells is above this cutoff")

# parameters for reference model: 
parser.add_argument("--labels_key_reference",
                    default=None,
                    help="")
parser.add_argument("--batch_key_reference",
                    default=None,
                    help="")
parser.add_argument("--layer_reference",
                    default=None,#or raw_counts
                    help="")
parser.add_argument("--categorical_covariate_keys_reference",
                    default=None,
                    help="")
parser.add_argument("--continuous_covariate_keys_reference",
                    default=None,
                    help="")
parser.add_argument("--max_epochs_reference",
                    default=None,
                    help="")

# parameters for spatial model: 
parser.add_argument("--labels_key_st",
                    default=None,
                    help="")
parser.add_argument("--batch_key_st",
                    default=None,
                    help="")
parser.add_argument("--layer_st",
                    default=None,#or raw_counts
                    help="")
parser.add_argument("--categorical_covariate_keys_st",
                    default=None,
                    help="")
parser.add_argument("--continuous_covariate_keys_st",
                    default=None,
                    help="")
parser.add_argument("--N_cells_per_location",
                    default=None,
                    help="")
parser.add_argument("--detection_alpha",
                    default=None,
                    help="")
parser.add_argument("--max_epochs_st",
                    default=None,
                    help="")  


args, opt = parser.parse_known_args()

L.info("running with args:")
L.debug(args)

figdir = args.figdir
if not os.path.exists(figdir):
    os.mkdir(figdir)
sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))

output_dir = args.output_dir
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


if args.N_cells_per_location is None: 
    L.error("The parameter 'N_cells_per_location' is not specified. The parameters 'N_cells_per_location' and 'detection_alpha' need to be specified!")
    sys.exit("The parameter 'N_cells_per_location' is not specified. The parameters 'N_cells_per_location' and 'detection_alpha' need to be specified!")
if args.detection_alpha is None: 
    L.error("The parameter 'detection_alpha' is not specified. The parameters 'N_cells_per_location' and 'detection_alpha' need to be specified!")
    sys.exit("The parameter 'detection_alpha' is not specified. The parameters 'N_cells_per_location' and 'detection_alpha' need to be specified!")
    
    
if (args.save_models is False) or (args.save_models == "False"): 
    save_models = False
else:
    save_models = True
      
if (args.remove_mt is True) or (args.remove_mt == "True"): 
    remove_mt = True
else:
    remove_mt = False


if args.categorical_covariate_keys_reference is not None: 
    categorical_covariate_keys_reference = list(args.categorical_covariate_keys_reference.split(","))
else:
    categorical_covariate_keys_reference = None
    
if args.categorical_covariate_keys_st is not None: 
    categorical_covariate_keys_st = list(args.categorical_covariate_keys_st.split(","))
else:
    categorical_covariate_keys_st = None
    
if args.continuous_covariate_keys_reference is not None: 
    continuous_covariate_keys_reference = list(args.continuous_covariate_keys_reference.split(","))
else:
    continuous_covariate_keys_reference = None
    
if args.continuous_covariate_keys_st is not None: 
    continuous_covariate_keys_st = list(args.continuous_covariate_keys_st.split(","))
else:
    continuous_covariate_keys_st = None

    
if args.max_epochs_reference is None: 
    max_epochs_reference = None
else: 
    max_epochs_reference = int(args.max_epochs_reference)
    
if args.max_epochs_st is None: 
    max_epochs_st = None
else: 
    max_epochs_st= int(args.max_epochs_st)




#1. read in the data
#spatial: 
mdata_spatial = mu.read(args.input_spatial)
adata_st = mdata_spatial.mod['spatial']
#single-cell: 
mdata_singlecell = mu.read(args.input_singlecell)
adata_sc = mdata_singlecell.mod['rna']



#2. Perform gene selection:
if args.gene_list != "None": # read in csv and subset both anndatas
    reduced_gene_set = pd.read_csv(args.gene_list, header = 0)
    reduced_gene_set.columns = ["HVGs"]
    adata_sc.var["selected_gene"] = adata_sc.var.index.isin(reduced_gene_set["HVGs"])
    adata_st.var["selected_gene"] = adata_st.var.index.isin(reduced_gene_set["HVGs"])
    adata_sc = adata_sc[:, adata_sc.var["selected_gene"]]
    adata_st = adata_st[:, adata_st.var["selected_gene"]]
    # check whether all genes are present in both, spatial & reference
    if set(adata_st.var.index) != set(adata_sc.var.index):
        L.error(
            "Not all genes of the gene list %s are present in the reference as well as in the ST data. Please provide a gene list where all genes are present in both, reference and ST.", args.gene_list)
        sys.exit(
            "Not all genes of the gene list are present in the reference as well as in the ST data. Please provide a gene list where all genes are present in both, reference and ST.")

else: # perform feature selection according to cell2loc
    if remove_mt is True: 
        adata_st.var["MT_gene"] = [gene.startswith("MT-") for gene in adata_st.var.index]
        adata_st.obsm["MT"] = adata_st[:, adata_st.var["MT_gene"].values].X.toarray()
        adata_st = adata_st[:, ~adata_st.var["MT_gene"].values]
    # intersect vars of reference and spatial
    shared_features = [feature for feature in adata_st.var_names if feature in adata_sc.var_names]
    adata_sc = adata_sc[:, shared_features]
    adata_st = adata_st[:, shared_features]
    # select features
    selected = cell2loc_filter_genes(adata_sc, figdir + "/gene_filter.png", cell_count_cutoff=float(args.cell_count_cutoff),
                                               cell_percentage_cutoff2=float(args.cell_percentage_cutoff2),
                                                nonz_mean_cutoff=float(args.nonz_mean_cutoff))

    adata_sc = adata_sc[:, selected]
    adata_st = adata_st[:, selected]

    

# 3. Fit regression model 
c2l.models.RegressionModel.setup_anndata(adata=adata_sc, 
                                         labels_key = args.labels_key_reference,
                                         layer= args.layer_reference, 
                                         batch_key= args.batch_key_reference,
                                         categorical_covariate_keys = categorical_covariate_keys_reference,
                                         continuous_covariate_keys =  continuous_covariate_keys_reference)
model_ref = c2l.models.RegressionModel(adata_sc)
model_ref.train(max_epochs=max_epochs_reference)

# plot elbo
cell2loc_plot_history(model_ref, figdir + "/ELBO_reference_model.png")

# export results
adata_sc = model_ref.export_posterior(adata_sc)
if "means_per_cluster_mu_fg" in adata_sc.varm.keys():
    inf_aver = adata_sc.varm["means_per_cluster_mu_fg"][[f"means_per_cluster_mu_fg_{i}" for i in adata_sc.uns["mod"]["factor_names"]]].copy()
else:
    inf_aver = adata_sc.var[[f"means_per_cluster_mu_fg_{i}" for i in adata_sc.uns["mod"]["factor_names"]]].copy()
inf_aver.columns = adata_sc.uns["mod"]["factor_names"]
inf_aver.to_csv(output_dir+"/Cell2Loc_inf_aver.csv")

# plot QC
cell2loc_plot_QC_reference(model_ref, figdir + "/QC_reference_reconstruction_accuracy.png", figdir + "/QC_reference_expression signatures_vs_avg_expression.png")

# save model and update mudata
mdata_singlecell.mod["rna"] = adata_sc
mdata_singlecell.update()
if save_models is True:
    model_ref.save(output_dir +"/Reference_model", overwrite=True)


       
# 4. Fit mapping model   
c2l.models.Cell2location.setup_anndata(adata=adata_st, 
                                         labels_key = args.labels_key_st,
                                         layer= args.layer_st, 
                                         batch_key= args.batch_key_st,
                                         categorical_covariate_keys = categorical_covariate_keys_st,
                                         continuous_covariate_keys =  continuous_covariate_keys_st)

        
model_spatial = c2l.models.Cell2location(adata = adata_st, cell_state_df=inf_aver, 
                                        N_cells_per_location=float(args.N_cells_per_location),
                                        detection_alpha=float(args.detection_alpha))
model_spatial.train(max_epochs=max_epochs_st)

#plot elbo
cell2loc_plot_history(model_spatial, figdir + "/ELBO_spatial_model.png")
#extract posterior
adata_st = model_spatial.export_posterior(adata_st)
#plot QC
cell2loc_plot_QC_reconstr(model_spatial, figdir + "/QC_spatial_reconstruction_accuracy.png")


#plot output
adata_st.obs[adata_st.uns["mod"]["factor_names"]] = adata_st.obsm["q05_cell_abundance_w_sf"]
sc.pl.spatial(adata_st,color=adata_st.uns["mod"]["factor_names"], show = False, save = "_Cell2Loc_q05_cell_abundance_w_sf.png") 


# save model and update mudata
mdata_spatial.mod["spatial"] = adata_st
mdata_spatial.update()
if save_models is True: 
    model_spatial.save(output_dir+"/Spatial_mapping_model", overwrite=True)


#6. save mudatas 
mdata_singlecell.write(output_dir+"/Cell2Loc_screference_output.h5mu")
mdata_spatial.write(output_dir+"/Cell2Loc_spatial_output.h5mu")


L.info("Done")

