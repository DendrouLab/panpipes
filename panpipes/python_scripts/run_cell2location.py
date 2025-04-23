'''
Run Cell2Location
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import cell2location as c2l
import scanpy as sc
import pandas as pd
import spatialdata as sd
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
parser.add_argument("--export_gene_by_spot",
                    default=False,
                    help="whether to save a gene by spot matrix for each cell type in a layer")


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
parser.add_argument("--accelerator_reference",
                    default="auto",
                    help="")
#parser.add_argument("--use_gpu_reference",
#                    default=True,
#                    help="")

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
parser.add_argument("--accelerator_spatial",
                    default="auto",
                    help="")
#parser.add_argument("--use_gpu_st",
#                    default=True,
#                    help="")


args, opt = parser.parse_known_args()

L.info("Running with params: %s", args)

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

if (args.export_gene_by_spot is False) or (args.export_gene_by_spot == "False"): 
    export_gene_by_spot = False
else:
    export_gene_by_spot = True
      
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


#if (args.use_gpu_reference is True) or (args.use_gpu_reference == "True"): 
#    use_gpu_reference = True
#else:
#    use_gpu_reference = False
#if (args.use_gpu_st is True) or (args.use_gpu_st == "True"): 
#    use_gpu_st = True
#else:
#    use_gpu_st = False



#1. read in the data
#spatial: 
L.info("Reading in SpatialData from '%s'" % args.input_spatial)
sdata_st = sd.read_zarr(args.input_spatial)
#mdata_spatial = mu.read(args.input_spatial)
#adata_st = mdata_spatial.mod['spatial']
#single-cell: 
L.info("Reading in reference MuData from '%s'" % args.input_singlecell)
mdata_singlecell = mu.read(args.input_singlecell)
adata_sc = mdata_singlecell.mod['rna']



#2. Perform gene selection:
if args.gene_list != "None": # read in csv and subset both anndatas
    if os.path.exists(args.gene_list):
        L.info("Reading in gene list file from '%s'" % args.gene_list)
    else: 
        L.error("File '%s' could not be found." % args.gene_list)
        sys.exit("File '%s' could not be found." % args.gene_list)
    reduced_gene_set = pd.read_csv(args.gene_list, header = 0)
    reduced_gene_set.columns = ["HVGs"]
    L.info("Subsetting data on gene list")
    adata_sc.var["selected_gene"] = adata_sc.var.index.isin(reduced_gene_set["HVGs"])
    sdata_st["table"].var["selected_gene"] = sdata_st["table"].var.index.isin(reduced_gene_set["HVGs"])
    adata_sc = adata_sc[:, adata_sc.var["selected_gene"]]
    sdata_st["table"] = sdata_st["table"][:, sdata_st["table"].var["selected_gene"]]
    # check whether all genes are present in both, spatial & reference
    if set(sdata_st["table"].var.index) != set(adata_sc.var.index):
        L.error(
            "Not all genes of the gene list %s are present in the reference as well as in the ST data. Please provide a gene list where all genes are present in both, reference and ST.", args.gene_list)
        sys.exit(
            "Not all genes of the gene list are present in the reference as well as in the ST data. Please provide a gene list where all genes are present in both, reference and ST.")

else: # perform feature selection according to cell2loc
    if remove_mt is True: 
        L.info("Removing MT genes")
        sdata_st["table"].var["MT_gene"] = [gene.startswith("MT-") for gene in sdata_st["table"].var.index]
        sdata_st["table"].obsm["MT"] = sdata_st["table"][:, sdata_st["table"].var["MT_gene"].values].X.toarray()
        sdata_st["table"] = sdata_st["table"][:, ~sdata_st["table"].var["MT_gene"].values]
    # intersect vars of reference and spatial
    L.info("Intersecting vars of reference and spatial ")
    shared_features = [feature for feature in sdata_st["table"].var_names if feature in adata_sc.var_names]
    adata_sc = adata_sc[:, shared_features]
    sdata_st["table"] = sdata_st["table"][:, shared_features]
    # select features
    L.info("Selecting features using 'cell2location.utils.filtering.filter_genes() function'")
    selected = cell2loc_filter_genes(adata_sc, figdir + "/gene_filter.png", cell_count_cutoff=float(args.cell_count_cutoff),
                                               cell_percentage_cutoff2=float(args.cell_percentage_cutoff2),
                                                nonz_mean_cutoff=float(args.nonz_mean_cutoff))
    L.info("Subsetting data on selected features")
    adata_sc = adata_sc[:, selected]
    sdata_st["table"] = sdata_st["table"][:, selected]

    

# 3. Fit regression model 
L.info("Setting up AnnData for the reference model")
c2l.models.RegressionModel.setup_anndata(adata=adata_sc, 
                                         labels_key = args.labels_key_reference,
                                         layer= args.layer_reference, 
                                         batch_key= args.batch_key_reference,
                                         categorical_covariate_keys = categorical_covariate_keys_reference,
                                         continuous_covariate_keys =  continuous_covariate_keys_reference)
model_ref = c2l.models.RegressionModel(adata_sc)
L.info("Training the reference model")
model_ref.train(max_epochs=max_epochs_reference, **{"accelerator": args.accelerator_reference})#use_gpu = use_gpu_reference)

# plot elbo
L.info("Plotting ELBO")
cell2loc_plot_history(model_ref, figdir + "/ELBO_reference_model.png")

# export results
L.info("Extracting the posterior of the reference model")
adata_sc = model_ref.export_posterior(adata_sc)
if "means_per_cluster_mu_fg" in adata_sc.varm.keys():
    inf_aver = adata_sc.varm["means_per_cluster_mu_fg"][[f"means_per_cluster_mu_fg_{i}" for i in adata_sc.uns["mod"]["factor_names"]]].copy()
else:
    inf_aver = adata_sc.var[[f"means_per_cluster_mu_fg_{i}" for i in adata_sc.uns["mod"]["factor_names"]]].copy()
inf_aver.columns = adata_sc.uns["mod"]["factor_names"]
inf_aver.to_csv(output_dir+"/Cell2Loc_inf_aver.csv")

# plot QC
L.info("Plotting QC plots")
cell2loc_plot_QC_reference(model_ref, figdir + "/QC_reference_reconstruction_accuracy.png", figdir + "/QC_reference_expression signatures_vs_avg_expression.png")

# save model 
if adata_sc.var.index.names[0] in adata_sc.var.columns: 
	adata_sc.var.index.names = [None]
mdata_singlecell.mod["rna"] = adata_sc
mdata_singlecell.update()
if save_models is True:
    L.info("Saving reference model to '%s'" % output_dir)
    model_ref.save(output_dir +"/Reference_model", overwrite=True)


       
# 4. Fit mapping model   
L.info("Setting up AnnData for the spatial model")
c2l.models.Cell2location.setup_anndata(adata=sdata_st["table"], 
                                         labels_key = args.labels_key_st,
                                         layer= args.layer_st, 
                                         batch_key= args.batch_key_st,
                                         categorical_covariate_keys = categorical_covariate_keys_st,
                                         continuous_covariate_keys =  continuous_covariate_keys_st)

        
model_spatial = c2l.models.Cell2location(adata = sdata_st["table"], cell_state_df=inf_aver, 
                                        N_cells_per_location=float(args.N_cells_per_location),
                                        detection_alpha=float(args.detection_alpha))
L.info("Training the spatial model")
model_spatial.train(max_epochs=max_epochs_st, **{"accelerator": args.accelerator_spatial})# use_gpu = use_gpu_st)

#plot elbo
L.info("Plotting ELBO")
cell2loc_plot_history(model_spatial, figdir + "/ELBO_spatial_model.png")
#extract posterior
L.info("Extracting the posterior of the spatial model")
sdata_st["table"] = model_spatial.export_posterior(sdata_st["table"])
#plot QC
L.info("Plotting QC plots")
cell2loc_plot_QC_reconstr(model_spatial, figdir + "/QC_spatial_reconstruction_accuracy.png")

# export a gene by spot matrix for each cell type
if export_gene_by_spot:
    # Compute expected expression per cell type
    expected_dict = model_spatial.module.model.compute_expected_per_cell_type(model_spatial.samples["post_sample_q05"], model_spatial.adata_manager)
    # Add to anndata layers
    for i, n in enumerate(model_spatial.factor_names_):
        sdata_st["table"].layers[n] = expected_dict['mu'][i]


#plot output
L.info("Plotting spatial embedding plot coloured by 'q05_cell_abundance_w_sf'")
sdata_st["table"].obs[sdata_st["table"].uns["mod"]["factor_names"]] = sdata_st["table"].obsm["q05_cell_abundance_w_sf"]
sc.pl.spatial(sdata_st["table"],color=sdata_st["table"].uns["mod"]["factor_names"], show = False, save = "_Cell2Loc_q05_cell_abundance_w_sf.png") 


# save model 
if sdata_st["table"].var.index.names[0] in sdata_st["table"].var.columns: 
	sdata_st["table"].var.index.names = [None]
#mdata_spatial.mod["spatial"] = adata_st
#mdata_spatial.update()
if save_models is True: 
    L.info("Saving spatial model to '%s'" % output_dir)
    model_spatial.save(output_dir+"/Spatial_mapping_model", overwrite=True)


#6. save data 
# check column names of spatial data
if any([' ' in col for col in sdata_st["table"].obs.columns]):
    L.warning("Spaces were found in column names of sdata['table']. Replacing spaces by underscores.")
    sdata_st["table"].obs.columns = sdata_st["table"].obs.columns.str.replace(' ', '_')
    
L.info("Saving SpatialData and MuData to '%s'" % output_dir)
mdata_singlecell.write(output_dir+"/Cell2Loc_screference_output.h5mu")
sdata_st.write(output_dir+"/Cell2Loc_spatial_output.zarr")


L.info("Done")

