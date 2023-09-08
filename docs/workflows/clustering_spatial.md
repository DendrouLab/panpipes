Clustering spatial data
=========================

Panpipes `clustering` workflow accepts both cell suspension datasets as well as spatial trascriptomics inputs that have been ingested with `qc_spatial` workflow and optionally filtered with the `spatial_preprocess` workflow.

`clustering` expects as a input a **single mudata object** with the `spatial` layer populated. 
Set `spatial: True` in the configuration file and customize the spatial modality clustering parameters exactly as you would for a single cell experiment. For more information check the default [clustering workflow](./clustering.md)
