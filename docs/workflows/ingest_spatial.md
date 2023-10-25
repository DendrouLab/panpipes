Ingesting spatial data
========================


Similar to the cell suspension workflow, Panpipes' `spatial_qc` ingests `Vizgen` and/or `Visium` data and saves the data into `MuData` objects. 
A main difference to the cell suspension `ingestion` workflow is that we are not concatenating the input data into a single matrix, but keeping the samples as separate `MuData` objects, each with a `spatial` layer. This is to ensure that the processing is not introducing any technical batch effect when tissue slides are very different in cell composition. In a future release, we will use [SpatialData](https://spatialdata.scverse.org/en/latest/tutorials/notebooks/notebooks.html) as a data format and framework to process multi-slides experiments.
 

### Steps

- Data is ingested into `MuData` objects with the modality `spatial`. The workflow generates one MuData per dataset.
    - Raw `MuData`s are saved into `./tmp`
- QC metrics are computed using `scanpy` functionalities: 
    - Basic QC metrics are computed using `sc.pp.calculate_qc_metrics`
    - (Optional) Compute cell-cycle scores using `sc.tl.score_genes_cell_cycle`. For that, the [default gene list](../../panpipes/resources/cell_cycle_genes.tsv) can be used or a path to a tsv file can be specified. 
    - (Optional) Custom genes actions. [Default gene list](../../panpipes/resources/qc_genelist_1.0.csv) can be used or a path to a csv file can be specified. 
        - Calculate proportions of gene groups, e.g. mitochondrial genes
        - Score genes using `sc.tl.score_genes`
    - `MuData`s with calculated QC metrics are saved in `qc.data`
    - Metadata (`.obs`) is saved into the current directory as tsv files
- Specified QC metrics are plotted in violin and spatial embedding plots
    - For `Vizgen` data, additonal histograms are plotted 




### Steps to run 

1.  Generate sample submission file. You can find more information about the generation [here](../usage/setup_for_spatial_workflows.md)
2.  (Optional) Generate QC gene lists as described in [Gene list format](../gene_list_format)
3.  Activate conda environment `conda activate pipeline_env`
4.  Generate yaml and log file (`panpipes spatial_qc config`)
5.  Edit the pipeline.yml file for your dataset
6.  Run complete QC pipeline with `panpipes spatial_qc make full --local`
7.  Use outputs to decide filtering thresholds
    -   **Note that the actual filtering occurs in the first step of the `preprocess` pipeline**


The [Ingesting 10X Visium data with Panpipes]() and [Ingesting MERFISH data with Panpipes]() tutorials guide you through the ingestion step by step. 


