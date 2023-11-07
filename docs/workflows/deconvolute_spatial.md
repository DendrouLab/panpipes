Deconvoluting spatial data
==========================

With the `deconvolution_spatial` workflow, one or multiple spatial slides can be deconvoluted in one run. For that, a `MuData` object for each slide is expected, with the spatial data saved in `mdata.mod["spatial"]`. The spatial slides are deconvoluted using the same reference. For the reference, one `MuData` with the gene expression data saved in `mdata.mod["rna"]` is expected as input. 

For now, the workflow provides the possibility to run deconvolution using `Cell2Location`.



## Steps
### Cell2Location 
For the reference and each spatial slide the following steps are run. **Note, that the same parameter setting is used for each slide.** 

- Gene selection. There are two possibilities for the gene selection: 
    - Genes of a user-provided feature set (csv-file) are used for deconvolution
    - Feature selection performed according to Cell2Location, i.e. via the function: `cell2location.utils.filtering.filter_genes`
    - **Note**: if no csv-file is provided by the user, the workflow will run the feature selection via the function: `cell2location.utils.filtering.filter_genes`. Thus, gene selection is **not optional**.  
- Regression/reference model is fitted and a plot of the training history as well as QC plots are saved in the `./figures/Cell2Location` directory. Additionally, a csv-file `Cell2Loc_inf_anver.csv` with the estimated expression of every gene in every cell type is saved in `./cell2location.output`.
- (Optional) Reference model is saved in `./cell2location.output`
- Spatial mapping model is fitted. Training history and QC plots are saved in the `./figures/Cell2Location` directory. Plots of the spatial embedding coloured by `q05_cell_abundance_w_sf` is also saved in `./figures/Cell2Location`.
- (Optional) Spatial mapping model is saved in `./cell2location.output`
- `MuData` objects of the spatial slide and the reference are saved in `./cell2location.output`. The `MuData` object of the spatial slide contains the estimated cell type abundances.


## Steps to run

1.  Activate conda environment `conda activate pipeline_env`
2.  Generate yaml and log file `panpipes deconvolution_spatial config`
3.  Specify the parameter setting in the pipeline.yml file 
4.  Run complete deconvolution workflow with `panpipes deconvolution_spatial make full --local`

The [Deconvoluting spatial data using Cell2Location]() tutorial guides you through deconvolution workflow of `Panpipes` step by step. 


