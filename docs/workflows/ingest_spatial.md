# Ingesting spatial data

Similarly to the cell suspension workflowm panpipes `spatial_qc` ingests the spatial data from `vizgen` or `visium` inputs.
Multiple files can be supplied providing a sample_annotation tsv file. The [usage](../usage/setup_for_spatial_workflows.md) section contains  more information on input file formatting.

### Steps

- Data is ingested into a `mudata` object with `spatial` layer. The worklfow generates one mudata per sample
- QC metrics are computed using `scanpy` and `squidpy` functionalities, such as distribution of total transcripts per cell, unique transcripts per cell, transcripts per FOV and the volume of the segmented cells. 
- Metrics are plotted as violin plots and histograms, and overlayed over the tissue slide.


One main difference with the cell suspension `ingestion` workflow is that we are not concatenating the input data into a single matrix, but keeping the samples as separate mudata objects, each with a `spatial` layer. This is to ensure that the processing is not introducing any technical batch effect when tissue slides are very different in cell composition. In a future release, we will use [SpatialData](https://spatialdata.scverse.org/en/latest/tutorials/notebooks/notebooks.html) as a data format and framework to process multi-slides experiments.





