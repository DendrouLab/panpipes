#!/usr/bin/python
import argparse
import yaml


parser = argparse.ArgumentParser()
#parser.add_argument('--mode_dictionary',default=mode_dictionary,type=json.loads)
parser.add_argument('--mode_dictionary',default="",
                    help="accepts a yaml speficication of modalities dictionary", 
                    type=yaml.safe_load)
parser.add_argument('--sample_id',
                    default=None,
                    help='')
parser.add_argument('--gex_infile',
                    default=None,
                    help='')
parser.add_argument('--gex_filetype',
                    default=None,
                    help='')
parser.add_argument('--adt_infile', 
                    default=None,
                    help='')
parser.add_argument('--adt_filetype', 
                    default=None,
                    help='')
parser.add_argument('--subset_adt_barcodes_to_gex', 
                    default=True,
                    help='')
parser.add_argument('--atac_infile', 
                    default=None,
                    help='')
parser.add_argument('--atac_filetype', 
                    default=None,
                    help='')
parser.add_argument('--tcr_filtered_contigs', 
                    default=None,
                    help='')
parser.add_argument('--tcr_filetype', 
                    default=None,
                    help='')
parser.add_argument('--bcr_filtered_contigs', 
                    default=None,
                    help='')
parser.add_argument('--bcr_filetype',
 default=None,
                    help='')
parser.add_argument('--output_file',
                    default=None,
                    help='')
parser.add_argument('--protein_var_table',
                    default=None,
                    help='')
parser.add_argument('--protein_new_index_col',
                    default=None,
                    help='')
parser.add_argument('--barcode_mtd_df',
                    default=None,
                    help='csv file contining barcode level metadata')
parser.add_argument('--barcode_mtd_metadatacols',
                    default=None,
                    help='comma separated strings listing the column you want to keep in barcode_mtd_df')                    
parser.add_argument('--per_barcode_metrics_file',
                    default=None,
                    help='ATAC/Multiome specific input file from csv')                    
parser.add_argument('--peak_annotation_file',
                    default=None,
                    help='ATAC/Multiome specific input file from csv')                    
parser.add_argument('--fragments_file',
                    default=None,
                    help='ATAC/Multiome specific input file from csv')                    


parser.set_defaults(verbose=True)
args, opt = parser.parse_known_args()

# import scanpy as sc
import pandas as pd
# import numpy as np
# from scipy.sparse import csr_matrix
import muon as mu
import warnings
from muon._atac.tools import add_peak_annotation, locate_fragments


import sys
import logging
L = logging.getLogger()
L.setLevel(logging.DEBUG)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)

from panpipes.funcs.io import check_filetype, load_mdata_from_multiple_files, update_intersecting_feature_names
from panpipes.funcs.processing import check_for_bool, add_var_mtd, update_var_index, intersection, intersect_obs_by_mod

L.debug(args)
# Make mdata object

# unimodal mu (check if all the modalities)
if isinstance(args.mode_dictionary, dict):
    mode_dictionary = args.mode_dictionary
else:
    mode_dictionary = yaml.safe_load(args.mode_dictionary) 

permf = [key for key, value in mode_dictionary.items() if value == True]
all_files = {"rna": [args.gex_infile, args.gex_filetype], 
            "prot": [args.adt_infile, args.adt_filetype], 
            "atac": [args.atac_infile, args.atac_filetype], 
            "tcr":[args.tcr_filtered_contigs, args.tcr_filetype],
            "bcr":[args.bcr_filtered_contigs, args.bcr_filetype]}
#subset to the modalities we want from permf
all_files = {nm: x  for (nm, x) in all_files.items() if nm in permf}
# this will work for any combination of modalities or number of modalities in all_files.
[check_filetype(x[0], x[1]) for x in all_files.values()]
mdata = load_mdata_from_multiple_files(all_files)
mdata.update()
mdata['rna'].var_names_make_unique()
mdata.update()
# now lets do some extra processing on the different modalities
if 'rna' in mdata.mod.keys():
    # keep only gex with at least 1 count (relevent when loading raw data)
    mdata['rna'].obs['total_counts'] = mdata['rna'].X.sum(axis=1) 
    mu.pp.filter_obs(mdata['rna'], 'total_counts', lambda x: (x > 1))
    mdata.update()
    # drop this columns as we don't need it anymore, and it will confuse things downstream
    mdata['rna'].obs = mdata['rna'].obs.drop(columns="total_counts")
if 'prot' in mdata.mod.keys():
    if 'rna' in mdata.mod.keys():
        if check_for_bool(args.subset_adt_barcodes_to_gex):
            intersect_obs_by_mod(mdata, ['rna', 'prot'])
    L.info(mdata['prot'].var.head())
    if args.protein_var_table is not None:
        try:
            df = pd.read_csv(args.protein_var_table, sep='\t', index_col=0)
            L.info("merging protein table with var")
            # add_var_mtd(mdata['prot'], df)
            mdata['prot'].var = mdata['prot'].var.merge(df, left_index=True, right_index=True)
            L.info(mdata['prot'].var.head())
            if args.protein_new_index_col is not None:
                L.info("updating prot.var index")
                update_var_index(mdata['prot'], args.protein_new_index_col)
            # we might want to split hashing antibodies into a separate modalities
            # we assume this has been inidicated in a "hashing_ab" column in the protein metadata file
            if "hashing_ab" in mdata['prot'].var.columns:
                # create new modality for hashing
                mdata.mod["hashing_ab"]=mdata["prot"][:, mdata["prot"].var["hashing_ab"]]
                # subset old modality to remove hashing
                mdata["prot"][:, ~mdata["prot"].var["hashing_ab"]]
        except FileNotFoundError:
            warnings.warn("protein metadata table not found")
if 'atac' in mdata.mod.keys():
    mdata['atac'].var_names_make_unique()
    if args.per_barcode_metrics_file is not None:
        L.info("Returning mudata with metadata from barcode metrics csv")
        per_barcode_metrics_file = pd.read_csv(args.per_barcode_metrics_file, index_col=0) #it's a csv by default and goes in atac obs
        mdata.mod['atac'].obs = mdata.mod['atac'].obs.merge(per_barcode_metrics_file, left_index=True, right_index=True)
    if args.peak_annotation_file is not None:
        L.info("Returning mudata with peak annotation in mdata['atac'].uns['atac']['peak_annotation']")
        add_peak_annotation(mdata, annotation=args.peak_annotation_file, sep="\t", return_annotation= False)
    if args.fragments_file is not None:
        L.info("Returning mudata with fragments file in mdata['atac'].uns['files']['fragments']")
        locate_fragments(mdata,fragments=args.fragments_file, return_fragments= False)

#make var names unique
for mm in mdata.mod.keys():
    mdata[mm].var_names_make_unique()

mdata.obs['sample_id'] = str(args.sample_id)

# copy the  to each modality
for mm in mdata.mod.keys():
    print("saving sample_id to each modality")
    # mdata[mm].obs['sample_id'] = mdata.obs['sample_id']
    mdata[mm].obs['sample_id'] = mdata.obs.loc[mdata[mm].obs_names,:]['sample_id']
    
mdata.update()

# # this is deliberatley in twice, for magic code reasons
# mdata.obs['sample_id'] = str(args.sample_id)

# make var
# add multiplexing data
if args.barcode_mtd_df is not None :
    L.info("Add barcode level metadata")
    # check demult_metadatacols exists and contains antibody
    barcode_metadatacols = args.barcode_mtd_metadatacols.split(',')
    # load the demultiplexing data
    barcode_mtd_df = pd.read_csv(args.barcode_mtd_df, index_col=0)
    mdata.obs = mdata.obs.merge(barcode_mtd_df, left_index=True, right_index=True)



L.info("saving data")
# this will write the file as anndata or muon depending on the class of mdata.
mdata.write(args.output_file)
L.info("done")
