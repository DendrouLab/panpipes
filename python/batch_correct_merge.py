import muon as mu
import argparse
import os

import sys
import logging
L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)


parser = argparse.ArgumentParser()
parser.add_argument('--preprocessed_mudata',
                    default='adata_scaled.h5mu',
                    help='')
parser.add_argument('--output_mudata',
                    default='mdata_corrected.h5mu',
                    help='')
parser.add_argument('--rna_correction_choice', default=None,
                    help='')
parser.add_argument('--prot_correction_choice', default=None,
                    help='')
parser.add_argument('--atac_correction_choice', default=None,
                    help='')
parser.add_argument('--multimodal_correction_choice', default=None,
                    help='')


args, opt = parser.parse_known_args()


L.info(args)

base_object = args.preprocessed_mudata

multi_mod = args.multimodal_correction_choice
if multi_mod is not None:
    base_object = "tmp/" + multi_mod + "_scaled_adata.h5mu"


    

uni_mod_paths = {}

if args.rna_correction_choice is not None:
    uni_mod_paths['rna']= "tmp/" + args.rna_correction_choice + "_scaled_adata_rna.h5ad"
if args.prot_correction_choice is not None:
    uni_mod_paths['prot']= "tmp/" + args.prot_correction_choice + "_scaled_adata_prot.h5ad"
if args.atac_correction_choice is not None:
    uni_mod_paths['atac'] = "tmp/" + args.atac_correction_choice + "_scaled_adata_atac.h5ad"
uni_mod_paths

base_obj = mu.read(base_object)


if multi_mod=="totalvi":
    totalvi_extra_obsm = {'prot': ['totalvi_protein_foreground_prob', 'totalvi_denoised_protein'], 'rna': [ 'totalvi_denoised_rna']}

    # we will want to keep these so extract them first
    totalvi_extra_obsm_keep = {}
    for k, v in totalvi_extra_obsm.items():
        print(k,v)
        totalvi_extra_obsm_keep[k] = {w:base_obj[k].obsm[w] for w in v}


# this overwrites the old modalities
for k,v in uni_mod_paths.items():
    base_obj.mod[k] = mu.read(v)

if multi_mod=="totalvi":
    # restore the totalvi extras
    for k, v in totalvi_extra_obsm_keep.items():
        for w in v:
            base_obj[k].obsm[ w] = totalvi_extra_obsm_keep[k][w]

base_obj.update()
base_obj.write(args.output_mudata)