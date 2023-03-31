import scanpy as sc
import glob
import pandas as pd
import argparse
import re
from  muon import MuData
import muon as mu

import warnings
# warnings.filterwarnings("ignore", category=UserWarning)
# warnings.filterwarnings("ignore", category=FutureWarning)

from panpipes.funcs.processing import  concat_mdatas

pd.set_option('display.max_rows', None)
# import os
# os.chdir("../../data/run_qc_general")
import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)



parser = argparse.ArgumentParser()
parser.add_argument('--input_files_str',
                    default='',
                    help='')
parser.add_argument('--output_file',
                    default='',
                    help='mdata.h5mu')
parser.add_argument("--submissionfile",
                    default="caf_example",
                    help="this file has all the samples of the experiment. it has at least 3 columns: \
                    sample_id,path and type. \
                    it can contain other anonymous columns which will be renamed using the metadatacols information")
parser.add_argument("--sampleprefix",
                    default="",
                    help="prefix to prepend when saving the metadata file")
parser.add_argument("--metadatacols",
                    default=None,
                    help="names of the metadata to rename the columns of the sample submission file ")
parser.add_argument("--join_type",
                    default="inner",
                    help="whether to concat anndata on the 'inner' (default) or 'outer'")
parser.add_argument('--barcode_mtd_df',
                    default=None,
                    help='csv file contining barcode level metadata')
parser.add_argument('--barcode_mtd_metadatacols',
                    default=None,
                    help='comma separated strings listing the column you want to keep in barcode_mtd_df')            
parser.add_argument('--protein_var_table',
                    default=None,
                    help='')
parser.add_argument('--protein_new_index_col',
                    default=None,
                    help='')        

parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()
L.info("running with:")
L.info(args)

if len(args.input_files_str) == 0:
    sys.exit("no input files detected")
    
lf = re.split(',', args.input_files_str)

sfile=args.submissionfile

#  lf = glob.glob("./tmp/*raw.h5mu")
# sfile= "CARTC_bonemarrow_samples.tsv"

caf = pd.read_csv(sfile, sep="\t")
# the modality argument is relevent only if use_muon is True.
# mdatas = [read_anndata(i, use_muon, modality="all") for i in lf]
mdatas = [mu.read(i) for i in lf]
# [x.var_names_make_unique() for x in mdatas]


if len(mdatas)==1:
    # if mdatas length is 0 then no concatenation required
    mdata = mdatas[0]
else:
    # @Fabiola I don't understand whats supposed to happen here
    # with the atac, should it just error out to say
    # that atac shouldn't be concatenated?
    temp = mdatas[0]
    if 'atac' in temp.mod.keys():
        mdata=mdatas[0]
        del temp
    elif 'prot' in temp.mod.keys() or 'rna' in temp.mod.keys():
        ## IF GEX and ADT is ok to concatenate ----------
        mdata = concat_mdatas(mdatas,
            batch_key="sample_id", 
            join_type=args.join_type)

# remove to avoid mem issues
del mdatas

### add metmdata from caf
L.debug(mdata.obs.columns)
L.debug(mdata.obs.head())

L.info("add metadata")
# add metmdata to each object
metadatacols = args.metadatacols
# metadatacols="mdata_cols"
# check sample_id is in metadatacols as a minimum requirement
if metadatacols == None :
    pass
elif metadatacols == 'None':
    pass
else:
    metadatacols = metadatacols.split(',')
    if 'sample_id' not in metadatacols:
        metadatacols = ['sample_id'] + metadatacols
    # check all metadatacols are in the caf
    if all([x in caf.columns for x in metadatacols]):
        # if mdata is actually a mudata it will be helpful 
        # (although annoying to store this in rna and protein as 
        # well since the batch columns etc will be here)
        # hopefullt his whole bit will become unnesacry because of 
        # https://github.com/scverse/mudata/issues/2
        for mod in mdata.mod.keys():
            L.debug(mod)
            L.debug(mdata[mod].obs.head())
            assay = mdata[mod]
            # make sure there is no index name because rep might have one as default
            assay.obs.index.name = None
            assay.obs = assay.obs.reset_index().merge(caf[metadatacols], on="sample_id", how="left").set_index('index').copy()
            # remove the index name introduced with the set index.
            assay.obs.index_name = None
            L.debug(assay.obs.head())
        mdata.update_obs()
        # merge metadata
        # obs = mdata.obs
        # obs['cellbarcode'] = obs.index  
        mdata.obs=mdata.obs.reset_index().merge(caf[metadatacols], on="sample_id", how="left").set_index('index')
        mdata.obs.index.name = None
        mdata.obs['cellbarcode'] = mdata.obs.index
        mdata.update_obs()
    else:
        # which column is missing?
        missingcols= " ".join([x for x in metadatacols if x not in caf.columns])
        L.error("Required columns missing form samples file: %s" % missingcols)
        sys.exit("Exiting because required columns not in samples file, check the log file for details")


# if 'rep' in mdata.mod.keys():
#     col_update = mdata['rep'].obs.columns[mdata['rep'].obs.columns.str.contains("count")]
#     mdata['rep'].obs[col_update] = mdata['rep'].obs[col_update].apply(pd.to_numeric)
#     mdata.update()
    
# add in cell level metadata
if args.barcode_mtd_df is not None :
    L.info("Add barcode level metadata")
    # check demult_metadatacols exists and contains antibody
    barcode_metadatacols = args.barcode_mtd_metadatacols.split(',')
    # load the demultiplexing data
    barcode_mtd_df = pd.read_csv(args.barcode_mtd_df, index_col=None)
    if len(barcode_mtd_df['sample_id'].unique().tolist()) > 1:
        barcode_mtd_df.index = barcode_mtd_df[['barcode_id', 'sample_id']].apply(lambda row: "-".join(row.values.astype(str)), axis=1)
    else:
        barcode_mtd_df.index = barcode_mtd_df['barcode_id']
    barcode_mtd_df= barcode_mtd_df.drop(columns=['sample_id'])
    # add for each modality
    for mod in mdata.mod.keys():
        mdata[mod].obs = mdata[mod].obs.merge(barcode_mtd_df, left_index=True, right_index=True)
        mdata.update()
    # add at top level
    mdata.obs = mdata.obs.merge(barcode_mtd_df, left_index=True, right_index=True)


L.debug(mdata.obs.columns)
L.debug(mdata.obs.head())
# update the protein variable to add in extra info like isotype and alternate name for adt
if args.protein_var_table is not None:
    try:
        df = pd.read_csv(args.protein_var_table, sep='\t', index_col=0)
        L.info("merging protein table with var")
        # add_var_mtd(mdata['prot'], df)
        var_df = mdata['prot'].var.merge(df, left_index=True, right_index=True)
        if args.protein_new_index_col is not None:
            L.info("updating prot.var index")
            # update_var_index(mdata['prot'], args.protein_new_index_col)
            var_df =var_df.reset_index().set_index(args.protein_new_index_col)
            var_df = var_df.rename(columns={'index':'orig_id'})
            var_df.index.name = None
        mdata['prot'].var = var_df
        mdata.update_var()
        mdata.update()
        # we might want to split hashing antibodies into a separate modalities
        # we assume this has been inidicated in a "hashing_ab" column in the protein metadata file
        if "hashing_ab" in mdata['prot'].var.columns:
            # create new modality for hashing
            mdata.mod["hashing_ab"]=mdata["prot"][:, mdata["prot"].var["hashing_ab"]]
            # subset old modality to remove hashing
            mdata.mod['prot'] = mdata["prot"][:, ~mdata["prot"].var["hashing_ab"]]
    except FileNotFoundError:
        warnings.warn("protein metadata table not found")
mdata.update()
# tidy up metadata
# move sample_id to the front
# cols = mdata.obs.columns.tolist()
# cols.insert(0, cols.pop(cols.index('sample_id')))
# mdata.obs = mdata.obs.reindex(columns=cols)
L.debug(mdata.obs.dtypes)

L.info("writing to file {}".format(str(args.output_file)))

mdata.write(args.output_file)

L.info("done")
