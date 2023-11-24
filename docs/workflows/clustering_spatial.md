Clustering spatial data
=========================

The `clustering` workflow accepts both cell suspension datasets as well as spatial transcriptomics data as input that have been ingested with the `qc_spatial` workflow and optionally filtered with the `spatial_preprocess` workflow.
The workflow expects a **single `MuData` object** with the spatial data saved in `mdata.mod["spatial"]`.
 
Set `spatial: True` in the configuration file and customize the spatial modality clustering parameters exactly as you would for a single cell experiment. For more information check the [clustering workflow](./clustering.md)


