QC
====



1. Generate sample submission file as described in :doc:`../setup_for_qc_mm` 

2. Generate qc genelists as described in :doc:`../gene_list_format`


3. For adt assay - generate the protein metadata file
   `example <(https://github.com/DendrouLab/panpipes/blob/main/resources/protein_metadata_w_iso.md)>`__.
   This file is integrated into the mdata[‘prot’].var slot.
4. Generate config file (``panpipes qc_mm config``)
5. Edit the pipeline.yml file for your dataset

   -  this is explained step by step within the pipeline.yml file

6. Run complete qc pipeline with ``panpipes qc_mm make full``
7. Use outputs to decide filtering thresholds.

   -  **Note that the actual filtering occurs in the first step of
      Preprocess pipeline**
   -  TODO: create doc to explain the pipeline outputs

The h5mu file outputted from ``qc_mm`` contains concatenated raw counts
from all samples in the submission file, plus qc metrics are computed,
and these qc metrics are visualised in a variety of plots to aid the
user to determine data quality and filtering thresholds.




.. Running pipeline modules from different entry points.
.. ''''''''''''''''''''''''''''''''''
..  :doc:`../different_entry_points`
