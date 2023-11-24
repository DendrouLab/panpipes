Ingesting spatial data
========================


Similar to the cell suspension workflow, `spatial_qc` ingests `Vizgen` and/or `Visium` data and saves the data into `MuData` objects. 
A primary difference to the cell suspension `ingestion` workflow is that we are not concatenating the input data into a single matrix, but keeping the samples as separate `MuData` objects, each with a `spatial` layer. This ensures that the processing does not introduce any technical batch effect when tissue slides are very different in cell composition. In a future release, we will use [SpatialData](https://spatialdata.scverse.org/en/latest/tutorials/notebooks/notebooks.html) as a data format and framework to process multi-slides experiments.
 

## Steps

- Data is ingested into `MuData` objects with the modality `spatial`. The workflow generates one MuData per dataset.
    - Raw `MuData` objects are saved into `./tmp`
- QC metrics are computed using `scanpy` functionalities: 
    - Basic QC metrics are computed using `sc.pp.calculate_qc_metrics`
    - (Optional) Compute cell-cycle scores using `sc.tl.score_genes_cell_cycle`. For that, the [default gene list](../../panpipes/resources/cell_cycle_genes.tsv) can be used or a path to a tsv file can be specified. 
    - (Optional) Custom genes actions. [Default gene list](../../panpipes/resources/qc_genelist_1.0.csv) can be used or a path to a csv file can be specified. 
        - Calculate proportions of gene groups, e.g. mitochondrial genes
        - Score genes using `sc.tl.score_genes`
    - `MuData` objects with calculated QC metrics are saved in `qc.data`
    - Metadata (`.obs`) is saved into the current directory as tsv files
- Specified QC metrics are plotted in violin and spatial embedding plots
    - For `Vizgen` data, additional histograms are plotted 




## Steps to run 

1.  Generate sample submission file. You can find more information about the generation [here](../usage/setup_for_spatial_workflows.md)
2.  (Optional) Generate QC gene lists as described in [gene list format](../usage/gene_list_format.md)
3.  Activate conda environment `conda activate pipeline_env`
4.  Generate yaml and log file `panpipes qc_spatial config`
5.  Specify the parameter setting in the pipeline.yml file
6.  Run complete QC workflow with `panpipes qc_spatial make full --local`
7.  Use outputs to decide filtering thresholds
    -   **Note that the actual filtering occurs in the first step of the `preprocess_spatial` workflow**


The [Ingesting 10X Visium data with Panpipes](https://panpipes-tutorials.readthedocs.io/en/latest/ingesting_visium_data/Ingesting_visium_data_with_panpipes.html) and [Ingesting MERFISH data with Panpipes](https://panpipes-tutorials.readthedocs.io/en/latest/ingesting_merfish_data/Ingesting_merfish_data_with_panpipes.html) tutorials guide you through the ingestion step by step. 


