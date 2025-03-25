'''
Concatenate spatial transcriptomics data
'''
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


import scanpy as sc
import spatialdata as sd

import os
import argparse
import sys
import logging
import glob

L = logging.getLogger()
L.setLevel(logging.INFO)
log_handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s: %(levelname)s - %(message)s')
log_handler.setFormatter(formatter)
L.addHandler(log_handler)
L.debug("testing logger works")



sc.settings.verbosity = 3

parser = argparse.ArgumentParser()

parser.add_argument("--input_dir",
                    default="./tmp/",
                    help="")
parser.add_argument("--output_dir",
                    default="./concatenated.data/",
                    help="")




args, opt = parser.parse_known_args()

L.info("Running with params: %s", args)

if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

L.info("Reading in all SpatialDatas from '%s'" % args.input_dir)
sdatas = []
for file in glob.glob(glob.escape(args.input_dir) + "/*.zarr"):
    sdatas.append(sd.read_zarr(file))

sdata = sd.concatenate(sdatas, concatenate_tables=True)


L.info("Saving concatenated SpatialData to '%s'" % args.output_dir)
sdata.write(args.output_dir + "concatenated.zarr")

L.info("Done")

