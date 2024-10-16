import argparse
import os
import numpy as np
import scanpy as sc
import muon as mu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import panpipes.funcs as pnp

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--filtered_mudata",
                    default=None,
                    help="input mudata after rna filtering(if rna in modality dictionary)")  
parser.add_argument("--bg_mudata",
                    default=None,
                    help="the mdata containing the raw measurements")
parser.add_argument("--channel_col",
                    default=None,
                    help="")
parser.add_argument("--normalisation_methods", default="clr,dsb",
                   help="comma separated list of normalisation methods")
parser.add_argument("--clr_margin",
                    default=0,
                    help="run clr with margin calculated on cells (columns) or on features (rows)")    
parser.add_argument("--quantile_clipping",
                    default=False,
                    help="")  
parser.add_argument("--store_as_x",
                    default=None,
                    help="which normalization to store in the X slot. Default to dsb if specified in normalization method")   
parser.add_argument("--figpath",
                    default="./prot_figures",
                    help="where to save the files")  
parser.add_argument("--save_mtx",
                    default=False,
                    help="if per channel normalization is applied, where to save the mtx file")  
parser.add_argument("--save_mudata_path",
                    default=None,
                    help="")   
parser.add_argument("--run_pca",
                    default=False,
                    help="whether to run pca on protein modality")   
parser.add_argument("--n_pcs",
                    default=None,
                    help="number of components to calculate")
parser.add_argument("--pca_solver",
                    default="arpack",
                    help="which PCA solver to use")
parser.add_argument("--color_by",
                    default="sample_id",
                    help="which columns to fetch from the protein .obs slot")


args, opt = parser.parse_known_args()
save_mtx=pnp.pp.check_for_bool(args.save_mtx)

norm_methods = args.normalisation_methods.split(',')


L.info("Running with params: %s", args)

figdir = args.figpath
if not os.path.exists(figdir):
    os.mkdir(figdir)

sc.settings.figdir = figdir
sc.set_figure_params(scanpy=True, fontsize=14, dpi=300, facecolor='white', figsize=(5,5))


# load filtered data - if use_muon is True this will return a mudata object, else returns an mudata object
L.info("Reading in MuData from '%s'" % args.filtered_mudata)
try:
    all_mdata = mu.read(args.filtered_mudata)
except FileNotFoundError:
    L.error("filtered_mudata file not found")
    sys.exit("filtered_mudata file not found")

# load bg data
if 'dsb' in norm_methods:
    if args.bg_mudata is not None:
        L.info("Reading in background MuData from " + args.bg_mudata)
        try:
            all_mdata_bg = mu.read(args.bg_mudata)
            # subset to just the channel 
        except FileNotFoundError:
            L.error("Background MuData object not found in " + args.bg_mudata)
            sys.exit("Background MuData object not found in " + args.bg_mudata)
    else:
        L.error("You must specify a background MuData to run dsb, containing both rna and prot")
        sys.exit("You must specify a background MuData to run dsb, containing both rna and prot")

    # RNA and PROT layer must be present in the object to run dsb
    if 'rna' not in list(all_mdata_bg.mod.keys()):
      L.error("The MuData you specified does not contain rna modality. Cannot run dsb")
      sys.exit("The MuData you specified does not contain rna modality. Cannot run dsb")

    # the RNA and PROT must intersect - they won't if experimental design included some cells without PROT data
    L.info("Checking that only cells with PROT signal present are in the background object")
    if all_mdata_bg['rna'].shape[0] != all_mdata_bg['prot'].shape[0]:
        mu.pp.intersect_obs(mdata_bg)
        all_mdata_bg.update()
    
    # checking that the same proteins are in foreground and background (since foreground might have been filtered)
    L.info("Checking that the same proteins are in foreground and background")
    if len(all_mdata['prot'].var_names) != len(all_mdata_bg['prot'].var_names):
        mu.pp.filter_var(all_mdata_bg['prot'], all_mdata['prot'].var_names)
        all_mdata_bg.update()


# find isotypes columns 
L.info("Looking for isotypes column in .var['isotype']")
if "isotype" in all_mdata["prot"].var.columns:
    isotypes = list(all_mdata["prot"].var_names[all_mdata["prot"].var.isotype])
    # exclude isotypes from downstream analysis by setting them as the only not highly variable genes.
    L.info("Excluding isotypes from downstream analysis by setting them as the only not HVGs")
    all_mdata['prot'].var['highly_variable'] = ~all_mdata['prot'].var['isotype'] 
else:
    L.info("No isotype column found in .var")
    isotypes = None
    
L.info("The following isotypes were found:")
L.info(isotypes)


# save raw counts
L.info("Checking if raw data is available")
if pnp.scmethods.X_is_raw(all_mdata["prot"]):
    L.info("Saving raw counts from .X to .layers['raw_counts']")
    all_mdata["prot"].layers["raw_counts"] = all_mdata["prot"].X.copy()
elif "raw_counts" in all_mdata["prot"].layers :
    L.info(".layers['raw_counts'] already exists and copying it to .X")
    all_mdata["prot"].X = all_mdata["prot"].layers['raw_counts'].copy()
else:
    L.error("X is not raw data and raw_counts layer not found")
    sys.exit("X is not raw data and raw_counts layer not found")



if args.channel_col is not None:
    ## running per channel (parallelise in future)
    channel_cols = list(pnp.processing.mu_get_obs( all_mdata, 
        features=[args.channel_col],
        modalities=["prot"]).iloc[:,0].unique())
    #this line works only if the obs of the entire metadata is the same as the prot obs. which may not always be the case
    if not args.channel_col in all_mdata.obs.columns:
        print("adding channel col to mdata obs since it's only in the prot obs")
        all_mdata.obs.loc[all_mdata["prot"].obs.index,args.channel_col] = all_mdata["prot"].obs[args.channel_col]
    mdata_per_sample = {si: all_mdata[all_mdata.obs[args.channel_col] == si] for si in channel_cols}
    # this data is only necessary if running dsb
    if 'dsb' in norm_methods:
        mdata_bg_per_sample = {si: all_mdata_bg[all_mdata_bg.obs[args.channel_col] == si] for si in channel_cols}
    for si in channel_cols:
        # si = list(mdata_per_sample.keys())[0]
        L.info(si)
        mdata=mdata_per_sample[si]
        plot_features = list(mdata["prot"].var_names)
        plot_features.sort()
        plot_features = [x for x in plot_features if x not in isotypes]
        # first run clr
        if 'clr' in norm_methods:
            # make sure to start from raw counts
            mdata["prot"].X = mdata["prot"].layers["raw_counts"].copy()
            pnp.scmethods.run_prot_normalise(mdata=mdata, mdata_bg=None,
                    method="clr",
                    clr_margin=int(args.clr_margin))
            L.info("saving ridgeplot")
            pnp.plotting.ridgeplot(mdata["prot"], features=plot_features, layer="clr",  splitplot=6)
            plt.savefig(os.path.join(args.figpath, str(si) + "_clr_ridgeplot.png"))
            if isotypes is not None:
                pnp.plotting.ridgeplot(mdata["prot"], features=isotypes, layer="clr",  splitplot=1)
                plt.savefig(os.path.join(args.figpath, str(si) + "_clr_ridgeplot_isotypes.png"))
            # save out the data in what format?
            if save_mtx:
                pnp.io.write_10x_counts(mdata["prot"], os.path.join("prot_clr" , str(si) ), layer="raw_counts")
        # then run dsb
        if 'dsb' in norm_methods:
            mdata_bg = mdata_bg_per_sample[si]
            mdata = mdata_per_sample[si]
            # make sure to start from raw counts
            mdata["prot"].X = mdata["prot"].layers["raw_counts"].copy()
            # plot hiustogram of bg
            mdata_bg["rna"].obs["log10umi"] = np.array(np.log10(mdata_bg["rna"].X.sum(axis=1) + 1)).reshape(-1)
            # mu.pl.histogram(mdata_bg["rna"], ["log10umi"], bins=50)
            # plt.savefig(os.path.join(args.figpath, si + "_log10umi.png"))
            # run dsb
            # this stores a layer named after the method as well as overwriting X
            pnp.scmethods.run_prot_normalise(mdata=mdata, 
                                            mdata_bg=mdata_bg,
                                            method="dsb", 
                                            isotypes=isotypes)
            pnp.plotting.ridgeplot(mdata["prot"], features=plot_features, layer="dsb",  splitplot=6)
            plt.savefig(os.path.join(args.figpath, str(si) + "_dsb_ridgeplot.png"))
            if isotypes is not None:
                pnp.plotting.ridgeplot(mdata["prot"], features=isotypes, layer="dsb",  splitplot=1)
                plt.savefig(os.path.join(args.figpath, str(si) + "_dsb_ridgeplot_isotypes.png"))
            if save_mtx:
                pnp.io.write_10x_counts(mdata["prot"], os.path.join("prot_dsb" , str(si)), layer="raw_counts")
else:
    # run on all the data (not on a channel basis)
    # first run clr
     # do the ridgeplots for clr
    plot_features = list(all_mdata["prot"].var_names)
    plot_features.sort()
    plot_features = [x for x in plot_features if x not in isotypes]
    if 'clr' in norm_methods:
        # make sure to start from raw counts
        all_mdata["prot"].X = all_mdata["prot"].layers["raw_counts"].copy()
        # run normalise
        # this stores a layer named after the method as well as overwriting X
        L.info("Normalizing data with CLR")
        pnp.scmethods.run_prot_normalise(mdata=all_mdata, 
                mdata_bg=None,
                method="clr",
                clr_margin=int(args.clr_margin))
        L.info("Plotting Ridgeplot")
        pnp.plotting.ridgeplot(all_mdata["prot"], features=plot_features, layer="clr",  splitplot=6)
        plt.savefig(os.path.join(args.figpath, "clr_ridgeplot.png"))
        if isotypes is not None:
            pnp.plotting.ridgeplot(all_mdata["prot"], features=isotypes, layer="clr",  splitplot=1)
            plt.savefig(os.path.join(args.figpath, "clr_ridgeplot_isotypes.png"))
        # save out the data in what format?
        if save_mtx:
            pnp.io.write_10x_counts(all_mdata["prot"], os.path.join("prot_clr"), layer="clr")
    if 'dsb' in norm_methods:
        # make sure to start from raw counts
        all_mdata["prot"].X = all_mdata["prot"].layers["raw_counts"].copy()
        # then run dsb
        # comput log 10 umi
        all_mdata_bg["rna"].obs["log10umi"] = np.array(np.log10(all_mdata_bg["rna"].X.sum(axis=1) + 1)).reshape(-1)
        mu.pl.histogram(all_mdata_bg["rna"], ["log10umi"], bins=50)
        plt.savefig(os.path.join(args.figpath, "all_log10umi.png"))
        # this stores a layer named after the method as well as overwriting X
        L.info("Normalizing data with dsb")
        pnp.scmethods.run_prot_normalise(mdata=all_mdata, 
                mdata_bg=all_mdata_bg, 
                method="dsb",
                isotypes=isotypes) 
        # apply quantile clipping # discussed in FAQs https://cran.r-project.org/web/packages/dsb/vignettes/dsb_normalizing_CITEseq_data.html
        if args.quantile_clipping:
            L.info("Quantile clipping")
            all_mdata['prot'].layers['dsb'] = pnp.scmethods.quantile_clipping(all_mdata['prot'], layer="dsb", inplace=False)
            all_mdata['prot'].X = all_mdata['prot'].layers['dsb'].copy()
        # plot ridgeplot
        L.info("Plotting Ridgeplot")
        pnp.plotting.ridgeplot(all_mdata["prot"], features=plot_features, layer="dsb",  splitplot=6)
        plt.savefig(os.path.join(args.figpath, "dsb_ridgeplot.png"))
        if isotypes is not None:
            pnp.plotting.ridgeplot(all_mdata["prot"], features=isotypes, layer="dsb",  splitplot=1)
            plt.savefig(os.path.join(args.figpath, "dsb_ridgeplot_isotypes.png"))
        # save out the data in what format?
        if save_mtx:
            pnp.io.write_10x_counts(all_mdata["prot"], os.path.join("prot_dsb"), layer="dsb")
    
    if args.store_as_x is not None:
        all_mdata["prot"].X = all_mdata["prot"].layers[args.store_as_x]
        
    # run pca on X 
    # this basically makes no sense if you have a small panel of antibodies.
    if args.run_pca:
        L.info("Running PCA")

        if all_mdata['prot'].var.shape[0] < int(args.n_pcs):
            L.warning("You have less features (%s) than number of PCs (%t) you intend to calculate." % (all_mdata['prot'].var.shape[0], args.n_pcs))
            n_pcs = all_mdata['prot'].var.shape[0] - 1
            L.info("Setting number of PCS to %i" % int(n_pcs))
        else:
            n_pcs = int(args.n_pcs)
        sc.tl.pca(all_mdata['prot'], n_comps=n_pcs, 
                        svd_solver=args.pca_solver, 
                        random_state=0) 

        # do some plots!
        L.info("Plotting PCA")
        sc.pl.pca_variance_ratio(all_mdata['prot'], log=True, n_pcs=n_pcs, save=".png")
        
        
        col_variables = args.color_by.split(",")
        
        col_variables = [a.strip() for a in col_variables]
        
        col_use = [var for var in col_variables if var in all_mdata['prot'].obs.columns]

        sc.pl.pca(all_mdata['prot'], color=col_use, save = "_vars.png")
        sc.pl.pca_loadings(all_mdata['prot'], components="1,2,3,4,5,6", save = ".png")
        sc.pl.pca_overview(all_mdata['prot'], save = ".png")

    if args.save_mudata_path is not None:
        all_mdata.update()
        L.info("Saving updated MuData to '%s'" % args.save_mudata_path)
        all_mdata.write(args.save_mudata_path)


L.info("Done")
