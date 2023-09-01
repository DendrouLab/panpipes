import pytest
import pandas as pd
import types
# from scpipelines.funcs.processing import extract_parameter_from_fname
import panpipes.funcs as pnp
from panpipes.funcs.io import gen_load_anndata_jobs,update_cellranger_col
from panpipes.pipeline_ingest import  unfilt_file, raw_file
import os
from cgatcore import pipeline as P


PARAMS = P.get_parameters("../scpipelines/pipeline_ingest/pipeline.yml")
PARAMS['submission_file'] = os.path.dirname(__file__) + "/sample_caf.tsv"
PARAMS['sample_prefix'] = "test"

#test generators
def test_gen_load_anndata_jobs():
    caf = pd.read_csv(PARAMS['submission_file'], sep='\t')
    assert isinstance(gen_load_anndata_jobs(caf=caf), types.GeneratorType)
    assert  next(gen_load_anndata_jobs(caf=caf, mode_dictionary={"RNA":True, "ADT":True})) == ('../sample_1_rna/filtered_feature_bc_matrix',
                                                     './tmp/sample1.h5mu', 
                                                     'sample1', 
                                                     'cellranger', 
                                                     '../sample_1_prot/raw_feature_bc_matrix', 'cellranger', 
                                                     None,None, 
                                                     None,None,
                                                     None,None,
                                                     None) 
# all the nones are for repertoire and atac paths 

# test defaults
# this basically checks that the default yml has not changed
def test_filenames():
    assert unfilt_file() == "test_unfilt.h5mu"
    assert raw_file() == "test_raw.h5mu"
    assert update_cellranger_col("path", raw=False) == ("path/filtered_feature_bc_matrix", "cellranger")
    assert update_cellranger_col("path", raw=True) == ("path/raw_feature_bc_matrix", "cellranger")
    # test with temp file
    with open('raw_feature_bc_matrix.h5', 'w'): pass
    assert update_cellranger_col("", raw=True) == ("raw_feature_bc_matrix.h5", "10X_h5")
    os.remove('raw_feature_bc_matrix.h5')

