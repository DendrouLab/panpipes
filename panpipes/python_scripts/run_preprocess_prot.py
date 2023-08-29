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
                    help="")  
parser.add_argument("--bg_mudata",
                    default=None,
                    help="")
parser.add_argument("--channel_col",
                    default=None,
                    help="")
parser.add_argument("--normalisation_methods", default="clr,dsb",
                   help="comma separated list of normalisation methods")
parser.add_argument("--clr_margin",
                    default=0,
                    help="")    
parser.add_argument("--quantile_clipping",
                    default=False,
                    help="")  
parser.add_argument("--store_as_x",
                    default=None,
                    help="")   
parser.add_argument("--figpath",
                    default="./adt_figures",
                    help="")  
parser.add_argument("--save_mtx",
                    default=False,
                    help="")  
parser.add_argument("--save_mudata_path",
                    default=None,
                    help="")   
              

args, opt = parser.parse_known_args()
# args = argparse.Namespace(filtered_mudata='test.h5mu', bg_mudata='/well/cartography/users/zsj686/non_cart_projects/005-multimodal_scpipelines/ingest/test_raw.h5mu', channel_col=None, normalisation_methods='clr,dsb', clr_margin='0', quantile_clipping='True', figpath='./figures/adt', save_mtx=False, save_mudata_path='test.h5mu')
save_mtx=pnp.pp.check_for_bool(args.save_mtx)

norm_methods = args.normalisation_methods.split(',')


L.info(args)

# load filtered data - if use_umon is True this will return a mudata object, else returns an mudata object

L.info("reading filtered mudata object")
try:
    all_mdata = mu.read(args.filtered_mudata)
except FileNotFoundError:
    sys.exit("filtered_mudata file not found")

# load bg data
if 'dsb' in norm_methods:
    if args.bg_mudata is not None:
        L.info("reading bg file object")
        try:
            all_mdata_bg = mu.read(args.bg_mudata)
            # subset to just the channel 
        except FileNotFoundError:
            sys.exit("bg_mudata object not found")
    else:
        sys.exit("must specify a bg mudata to run dsb, containing both rna and prot")
    
    # checking that the samne proteins are in foreground and background (since foreground might have been filtered)
    if len(all_mdata['prot'].var_names) != len(all_mdata_bg['prot'].var_names):
        mu.pp.filter_var(all_mdata_bg['prot'], all_mdata['prot'].var_names)
        all_mdata_bg.update()


# find isotypes columns 
if "isotype" in all_mdata["prot"].var.columns:
    isotypes = list(all_mdata["prot"].var_names[all_mdata["prot"].var.isotype])
    # exclude isotypes from downstream analysis by setting them as the only not highly variable genes.
    all_mdata['prot'].var['highly_variable'] = ~all_mdata['prot'].var['isotype'] 
else:
    isotypes = None
    
L.info("isotypes found:")
L.info(isotypes)


# save raw counts
if pnp.scmethods.X_is_raw(all_mdata["prot"]):
    all_mdata["prot"].layers["raw_counts"] = all_mdata["prot"].X.copy()
elif "raw_counts" in all_mdata["prot"].layers :
    L.info("raw_counts layer already exists")
    all_mdata.X = all_mdata["prot"].layers['raw_counts'].copy()
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
            pnp.scmethods.run_adt_normalise(mdata=mdata, mdata_bg=None,
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
                pnp.io.write_10x_counts(mdata["prot"], os.path.join("adt_clr" , str(si) ), layer="raw_counts")
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
            pnp.scmethods.run_adt_normalise(mdata=mdata, 
                                            mdata_bg=mdata_bg,
                                            method="dsb", 
                                            isotypes=isotypes)
            pnp.plotting.ridgeplot(mdata["prot"], features=plot_features, layer="dsb",  splitplot=6)
            plt.savefig(os.path.join(args.figpath, str(si) + "_dsb_ridgeplot.png"))
            if isotypes is not None:
                pnp.plotting.ridgeplot(mdata["prot"], features=isotypes, layer="dsb",  splitplot=1)
                plt.savefig(os.path.join(args.figpath, str(si) + "_dsb_ridgeplot_isotypes.png"))
            if save_mtx:
                pnp.io.write_10x_counts(mdata["prot"], os.path.join("adt_dsb" , str(si)), layer="raw_counts")
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
        pnp.scmethods.run_adt_normalise(mdata=all_mdata, 
                mdata_bg=None,
                method="clr",
                clr_margin=int(args.clr_margin))
        pnp.plotting.ridgeplot(all_mdata["prot"], features=plot_features, layer="clr",  splitplot=6)
        plt.savefig(os.path.join(args.figpath, "clr_ridgeplot.png"))
        if isotypes is not None:
            pnp.plotting.ridgeplot(all_mdata["prot"], features=isotypes, layer="clr",  splitplot=1)
            plt.savefig(os.path.join(args.figpath, "clr_ridgeplot_isotypes.png"))
        # save out the data in what format?
        if save_mtx:
            pnp.io.write_10x_counts(all_mdata["prot"], os.path.join("adt_clr"), layer="clr")
    if 'dsb' in norm_methods:
        # make sure to start from raw counts
        all_mdata["prot"].X = all_mdata["prot"].layers["raw_counts"].copy()
        # then run dsb
        # comput log 10 umi
        all_mdata_bg["rna"].obs["log10umi"] = np.array(np.log10(all_mdata_bg["rna"].X.sum(axis=1) + 1)).reshape(-1)
        mu.pl.histogram(all_mdata_bg["rna"], ["log10umi"], bins=50)
        plt.savefig(os.path.join(args.figpath, "all_log10umi.png"))
        # this stores a layer named after the method as well as overwriting X
        pnp.scmethods.run_adt_normalise(mdata=all_mdata, 
                mdata_bg=all_mdata_bg, 
                method="dsb",
                isotypes=isotypes) 
        # apply quantile clipping # discussed in FAQs https://cran.r-project.org/web/packages/dsb/vignettes/dsb_normalizing_CITEseq_data.html
        if args.quantile_clipping:
            all_mdata['prot'].layers['dsb'] = pnp.scmethods.quantile_clipping(all_mdata['prot'], layer="dsb", inplace=False)
            all_mdata['prot'].X = all_mdata['prot'].layers['dsb'].copy()
        # plot ridgeplot
        pnp.plotting.ridgeplot(all_mdata["prot"], features=plot_features, layer="dsb",  splitplot=6)
        plt.savefig(os.path.join(args.figpath, "dsb_ridgeplot.png"))
        if isotypes is not None:
            pnp.plotting.ridgeplot(all_mdata["prot"], features=isotypes, layer="dsb",  splitplot=1)
            plt.savefig(os.path.join(args.figpath, "dsb_ridgeplot_isotypes.png"))
        # save out the data in what format?
        if save_mtx:
            pnp.io.write_10x_counts(all_mdata["prot"], os.path.join("adt_dsb"), layer="dsb")
    
    if args.store_as_x is not None:
        all_mdata["prot"].X = all_mdata["prot"].layers[args.store_as_x]
        
    # run pca on X 
    # this basically makes no sense if you have a small panel of antibodies.
    sc.tl.pca(all_mdata['prot'], 
              n_comps=min(50,all_mdata['prot'].var.shape[0]-1), 
              svd_solver='arpack', random_state=0) 

    if args.save_mudata_path is not None:
        all_mdata.update()
        all_mdata.write(args.save_mudata_path)


L.info("done")
