import argparse
import os
import numpy as np
import scanpy as sc
import muon as mu
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

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--filtered_mudata",
                    default=None,
                    help="")  
parser.add_argument("--raw_mudata",
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
save_mtx=pnp.pp.check_for_bool(args.save_mtx)

norm_methods = args.normalisation_methods.split(',')

L.info(args)
# load filtered data - if use_umon is True this will return a mudata object, else returns an mudata object

L.info("reading filtered mudata object")
try:
    all_mdata = mu.read(args.filtered_mudata)
except FileNotFoundError:
    sys.exit("filtered_mudata file not found")

# load raw data
if 'dsb' in norm_methods:
    if args.raw_mudata is not None:
        L.info("reading raw file object")
        try:
            all_mdata_raw = mu.read(args.raw_mudata)
            # subset to just the channel 
        except FileNotFoundError:
            sys.exit("raw_mudata object not found")
    else:
        sys.exit("must specify a raw mudata to run dsb, containing both rna and prot")


# find isotypes columns 
if "isotype" in all_mdata["prot"].var.columns:
    isotypes = list(all_mdata["prot"].var_names[all_mdata["prot"].var.isotype])
else:
    isotypes = None
    
L.info("isotypes found:")
L.info(isotypes)


# save raw counts
all_mdata["prot"].layers["raw_counts"] = all_mdata["prot"].X.copy()


if args.channel_col is not None:
    ## running per channel (parallelise in future)
    channel_cols = list(pnp.processing.mu_get_obs( all_mdata, 
        features=[args.channel_col],
        modalities=["prot"]).iloc[:,0].unique())
    mdata_per_sample = {si: all_mdata[all_mdata.obs[args.channel_col] == si] for si in channel_cols}
    # this data is only necessary if running dsb
    if 'dsb' in norm_methods:
        mdata_raw_per_sample = {si: all_mdata_raw[all_mdata_raw.obs[args.channel_col] == si] for si in channel_cols}
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
            mdata["prot"].X = mdata["prot"].layers["raw_counts"]
            pnp.scmethods.run_adt_normalise(mdata=mdata, mdata_raw=None,
                    method="clr",
                    clr_margin=int(args.clr_margin))
            L.info("saving ridgeplot")
            pnp.plotting.ridgeplot(mdata["prot"], features=plot_features, layer="clr",  splitplot=6)
            plt.savefig(os.path.join(args.figpath, si + "_clr_ridgeplot.png"))
            if isotypes is not None:
                pnp.plotting.ridgeplot(mdata["prot"], features=isotypes, layer="clr",  splitplot=1)
                plt.savefig(os.path.join(args.figpath, si + "_clr_ridgeplot_isotypes.png"))
            # save out the data in what format?
            if save_mtx:
                pnp.io.write_10x_counts(mdata["prot"], os.path.join("adt_clr" , si ), layer="raw_counts")
        # then run dsb
        if 'dsb' in norm_methods:
            mdata_raw = mdata_raw_per_sample[si]
            mdata = mdata_per_sample[si]
            # make sure to start from raw counts
            mdata["prot"].X = mdata["prot"].layers["raw_counts"]
            # plot hiustogram of bg
            mdata_raw["rna"].obs["log10umi"] = np.array(np.log10(mdata_raw["rna"].X.sum(axis=1) + 1)).reshape(-1)
            # mu.pl.histogram(mdata_raw["rna"], ["log10umi"], bins=50)
            # plt.savefig(os.path.join(args.figpath, si + "_log10umi.png"))
            # run dsb
            pnp.scmethods.run_adt_normalise(mdata=mdata, 
                                            mdata_raw=mdata_raw,
                                            method="dsb", 
                                            isotypes=isotypes)
            pnp.plotting.ridgeplot(mdata["prot"], features=plot_features, layer="dsb",  splitplot=6)
            plt.savefig(os.path.join(args.figpath, si + "_dsb_ridgeplot.png"))
            if isotypes is not None:
                pnp.plotting.ridgeplot(mdata["prot"], features=isotypes, layer="dsb",  splitplot=1)
                plt.savefig(os.path.join(args.figpath, si + "_dsb_ridgeplot_isotypes.png"))
            if save_mtx:
                pnp.io.write_10x_counts(mdata["prot"], os.path.join("adt_dsb" , si), layer="raw_counts")
else:
    # run on all the data (not on a channel basis)
    # first run clr
     # do the ridgeplots for clr
    plot_features = list(all_mdata["prot"].var_names)
    plot_features.sort()
    plot_features = [x for x in plot_features if x not in isotypes]
    if 'clr' in norm_methods:
        # make sure to start from raw counts
        all_mdata["prot"].X = all_mdata["prot"].layers["raw_counts"]
        # run normalise
        pnp.scmethods.run_adt_normalise(mdata=all_mdata, 
                mdata_raw=None,
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
        all_mdata["prot"].X = all_mdata["prot"].layers["raw_counts"]
        # then run dsb
        # comput log 10 umi
        all_mdata_raw["rna"].obs["log10umi"] = np.array(np.log10(all_mdata_raw["rna"].X.sum(axis=1) + 1)).reshape(-1)
        mu.pl.histogram(all_mdata_raw["rna"], ["log10umi"], bins=50)
        plt.savefig(os.path.join(args.figpath, "all_log10umi.png"))
        pnp.scmethods.run_adt_normalise(mdata=all_mdata, 
                mdata_raw=all_mdata_raw, 
                method="dsb",
                isotypes=isotypes) 
        # apply quantile clipping # discussed in FAQs https://cran.r-project.org/web/packages/dsb/vignettes/dsb_normalizing_CITEseq_data.html
        if args.quantile_clipping:
            all_mdata['prot'].layers['dsb_clipped'] = pnp.scmethods.quantile_clipping(all_mdata['prot'], layer="dsb", inplace=False)
            all_mdata['prot'].X = all_mdata['prot'].layers['dsb_clipped']
        # plot ridgeplot
        pnp.plotting.ridgeplot(all_mdata["prot"], features=plot_features, layer="dsb",  splitplot=6)
        plt.savefig(os.path.join(args.figpath, "dsb_ridgeplot.png"))
        if isotypes is not None:
            pnp.plotting.ridgeplot(all_mdata["prot"], features=isotypes, layer="dsb",  splitplot=1)
            plt.savefig(os.path.join(args.figpath, "dsb_ridgeplot_isotypes.png"))
        # save out the data in what format?
        if save_mtx:
            pnp.io.write_10x_counts(all_mdata["prot"], os.path.join("adt_dsb"), layer="dsb")

    # run pca on dsb normalised (if dsb was run otherwise clr is in X)
    sc.tl.pca(all_mdata['prot'], n_comps=50, svd_solver='arpack', random_state=0) 

    if args.save_mudata_path is not None:
        all_mdata.update()
        pnp.io.write_anndata(all_mdata, args.save_mudata_path, use_muon=True, modality='all')



L.info("Done")

