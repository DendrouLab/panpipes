<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style>

# Ingestion YAML

In this documentation, the parameters of the `ingest` configuration yaml file are explained. 
This file is generated running `panpipes ingest config`.  <br> The individual steps run by the pipeline are described in the [ingestion workflow](../workflows/qc.md). 

When running the ingestion workflow, panpipes provides a basic `pipeline.yml` file.
To run the workflow on your own data, you need to specify the parameters described below in the `pipeline.yml` file to meet the requirements of your data.
However, we do provide pre-filled versions of the `pipeline.yml` file for individual [tutorials](https://panpipes-pipelines.readthedocs.io/en/latest/tutorials/index.html).

For more information on functionalities implemented in `panpipes` to read the configuration files, such as reading blocks of parameters and reusing blocks with  `&anchors` and `*scalars`, please check [our documentation](./useful_info_on_yml.md)


You can download the different ingestion `pipeline.yml` files here:
- Basic `pipeline.yml` file (not prefilled) that is generated when calling `panpipes ingest config: [Download here](https://github.com/DendrouLab/panpipes/blob/main/panpipes/panpipes/pipeline_ingest/pipeline.yml)
- `pipeline.yml` file for [Ingesting data Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/ingesting_data/Ingesting_data_with_panpipes.html): [Download here](https://github.com/DendrouLab/panpipes-tutorials/blob/main/docs/ingesting_data/pipeline.yml)
- `pipeline.yml` file for [Ingesting Mouse data Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/ingesting_mouse/Ingesting_mouse_data_with_panpipes.html): [Download here](https://github.com/DendrouLab/panpipes-tutorials/blob/main/docs/ingesting_mouse/pipeline.yml)
- `pipeline.yml` file for [Ingesting multimodal (CITE-Seq + VDJ) data Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/ingesting_multimodal_data/ingesting_multimodal_data.html): [Download here](https://github.com/DendrouLab/panpipes-tutorials/blob/main/docs/ingesting_multimodal_data/pipeline.yml)
- `pipeline.yml` file for [Ingesting multiome data from cellranger Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/ingesting_multiome/ingesting_mome.html): [Download here](https://github.com/DendrouLab/panpipes-tutorials/blob/main/docs/ingesting_multiome/pipeline.yml)


## Compute resources options

<span class="parameter">resources</span><br>
Computing resources to use, specifically the number of threads used for parallel jobs, check [threads_tasks_panpipes.xlsx](https://github.com/DendrouLab/panpipes/blob/threads_doc_g/docs/yaml_docs/threads_tasks_panpipes.xlsx) for more information on which threads each specific task requires.
Specified by the following three parameters:
  - <span class="parameter">threads_high</span> `Integer`, Default: 1<br>
        Number of threads used for high intensity computing tasks. 
        For each thread, there must be enough memory to load all your input files at once and create the MuData object.

  - <span class="parameter">threads_medium</span> `Integer`, Default: 1<br>
        Number of threads used for medium intensity computing tasks.
        For each thread, there must be enough memory to load your mudata and do computationally light tasks.

  - <span class="parameter">threads_low</span> `Integer`, Default: 1<br>
  	    Number of threads used for low intensity computing tasks.
        For each thread, there must be enough memory to load text files and do plotting, requires much less memory than the other two.

<span class="parameter">condaenv</span> `String` (Path)<br>
    Path to conda environment that should be used to run panpipes.
    Leave blank if running native or your cluster automatically inherits the login node environment

## Loading and merging data options
### Project name and data format

<span class="parameter">project</span> `String`, Default: "test"<br>
    Project name.

<span class="parameter">sample_prefix</span> `String`, Default: "test"<br>
    Prefix for sample names.

<span class="parameter">use_existing_h5mu</span> `Boolean` (True/False), Default: False<br>
    If you want to read an existing MuData object (.h5mu file), set to True. Move the object to the folder where you intend to run the ingestion workflow, and call it ${sample_prefix}_unfilt.h5mu where ${sample_prefix} = sample_prefix argument above.

<span class="parameter">submission_file</span> `String`, Mandatory parameter<br>
    Submission file name (e.g. sample_file_qc.txt).
    Ensure that the submission file must be in the right format.
    For ingest, the submission file must contain the follwing columns: sample_id  rna_path  rna_filetype  (prot_path  prot_filetype tcr_path  tcr_filetype etc.)
    An example submission file can be found at resources/sample_file_mm.txt

<span class="parameter">metadatacols</span> `String` (comma-separated)<br>
    Specify the metadata columns from the submission file that you want to include in for the downstream processing of the data.
    The specified columns will we included in the created AnnData object that stores all data. 
    Provide the columns as a comma-separated String, for example: batch,disease,sex

<span class="parameter">concat_join_type</span> `String`, Default: inner<br>
    Specifies which concat join is performed on the MuData objects.
    We recommended `inner`.
    See the [AnnData documentation](https://anndata.readthedocs.io/en/latest/concatenation.html#inner-and-outer-joins) for details.


### Modalities in the project

It is crucial to specify which modalities are present in the data and need to be processed by panpipes.
The Quality Control (QC) is run for each modality independently.

<span class="parameter">modalities:</span><br>
    Specify which modalities are included in the data by setting the respective modality to True.
    Leave empty (None) or False to signal this modality is not part of the experiment.
    The modalities are processed in the order of the following list:
  - <span class="parameter">rna</span> `Boolean`, Default: False<br>

  - <span class="parameter">prot</span> `Boolean`, Default: False<br>

  - <span class="parameter">bcr</span> `Boolean`, Default: False<br>

  - <span class="parameter">tcr</span> `Boolean`, Default: False<br>

  - <span class="parameter">atac</span> `Boolean`, Default: False<br>
    
### Integrating barcode level data

<span class="parameter">barcode_mtd</span><br>
    Optionally, you can ingest cell-barcode-level metadata, such as demultiplexing with hashtags, chemical tags or lipid tagging.
    In case you have cell level metadata such as results from a demultiplexing algorithm, you can incorporate it into the MuData object.
    The demultiplexing results should be stored in one csv file containing 2 columns, barcode_id, and sample_id.
    Note that the sample_id should match the sample_id column in the submission file.
  
  - <span class="parameter">include</span> `Boolean`, Default: False<br>
        Set to True if you want to include barcode level data.

  - <span class="parameter">path</span> `String`<br>
        Insert path to the csv file containing the demultiplexing results.
           Requires include to be set to True.

  - <span class="parameter">metadatacols</span> `String`<br>
        Provide the metadata columns in the file specified via path that should be incorporated in the MuData object.
    

### Loading Protein data - additional options
Optionally, you can provide additional information on the Protein modality by specifying the following parameters.

<span class="parameter">protein_metadata_table</span> `String`<br>
    File containing additional information on the proteins.
    By default, the ingestion workflow will use the first column of the cellranger features.tsv.gz.
    If you want to add additional information about the antibodies (e.g. whether they are hashing antibodies or isotypes), you need to provide an additional table containing these information.
    The first column of the table must be equivalent to the first column of cellrangers features.tsv.gz and it must be stored as a txt file.
    For example, you could create a table that include a column called isotype (containing True or False) to run QC isotypes, and a column called hashing_ab (containing True or False) if your dataset contains hashing antibodies which you wish to split into a separate modality.

<span class="parameter">index_col_choice</span> `String`<br>
    Column name from protein_metadata_table that should be included in the MuData object.
    If you want to update the MuData object storing all our data with a column from the table specified in protein_metadata_table, specify the column name here.
    Ensure that there are no overlaps with the gene symbols used for the rna modality, otherwise you might face issues down the line.
    If there are overlaps with the rna gene symbols, then you might have trouble down the line
    Consider adding e.g. a prefix to make the protein labels unique (e.g. "prot_CD4").

<span class="parameter">load_prot_from_raw</span> `Boolean`, Default: False<br>
    In case separate prot and rna 10X outputs are provided, the pipeline will load the filtered rna 10X outputs and the raw prot 10X counts.

<span class="parameter">subset_prot_barcodes_to_rna</span> `Boolean`, Default: False<br>
    If providing separate prot and rna 10X output (`load_prot_from_raw`: True), consider setting `subset_prot_barcodes_to_rna` to True.
    Set to True if you want to treat filtered rna barcodes as the "real" cells (which we assume one wants to do in most cases), and, additionally, the out prot assay barcodes should match the rna ones.
    If set to False, the full prot matrix will be loaded into the MuData object.

## Quality Control (QC) options
### Processing of 10X cellranger files

<span class="parameter">plot_10X_metrics</span> `Boolean`, Default: False<br>
    Specify if you want to plot 10x cellranger QC metrics.
    Set this parameter to False if you're not using cellranger folders as input.
    If starting from cellranger output, you can parse multiple metrics_summary.csv files and plot them together.
    The ingestion workflow will search for the metrics_summary file(s) in the "outs" folder and will stop if it doesn't find it.

### Doublet detection on RNA modality
For doublet detection, we use the [Scrublet tool](https://www.sciencedirect.com/science/article/pii/S2405471218304745?via%3Dihub) for detection of doublets in the data. Doublets are cells that contain two distinct cell barcodes, and are often a result of two cells being encapsulated in the same droplet during the library preparation step.

<span class="parameter">scr</span><br>
    Specify the following parameters for doublet detection using the tool Scrublet.

  - <span class="parameter">expected_doublet_rate</span> `Float`, Default: 0.06<br>
        Fraction of observations (cells) expected to be doublets.
        Note that results are not particularly sensitive to this parameter.

  - <span class="parameter">sim_doublet_ratio</span> `Float`, Default: 2<br>
        Number of doublets to simulate, relative to the number of total observations.
        Note that setting this too high drastically increases computationally complexity.
        The minimum value we tested is 0.5.

  - <span class="parameter">n_neighbours</span> `Integer`, Default: 20<br>
        Number of neighbors used to construct the k-nearest-neighbor (KNN) classifier for the observations and simulated doublets.
        The value (round(0.5*sqrt(n_cells))) usually works well.
    
  - <span class="parameter">min_counts</span> `Integer`, Default: 2<br>
        Used for gene filtering prior to PCA.
        Genes expressed less than `min_counts` times in `min_cells` cells (see next parameter) are excluded from the data for further processing.

  - <span class="parameter">min_cells</span> `Integer`, Default: 3<br>
        Used for gene filtering prior to PCA.
        Genes expressed less than `min_counts` (see previous parameter) times in `min_cells` cells are excluded from the data for further processing."
  
  - <span class="parameter">min_gene_variability_pctl</span> `Integer`, Default: 85<br>
        Used for gene filtering prior to PCA.
        Keeps only the highly variable genes in the top min_gene_variability_pctl percentile for downstream processing, as measured by the v-statistic [Klein et al., Cell 2015].
  
  - <span class="parameter">n_prin_comps</span> `Integer`, Default: 30<br>
        Number of principal components used to embed the observations in a lower-dimensional space.
        PCA is run prior to k-nearest-neighbor graph construction.
  
  - <span class="parameter">use_thr</span> `Boolean`, Default: True<br>
        Specify if you want to use a user-defined threshold in plots to distinguish between true and false doublets.
        The actual threshold is defined in the following parameter.
        If set to False, doublet detection will be run with the default threshold calculated by Scrublet.
        Note that this threshold applies to plots only, meaning no actual filtering takes place here.
  
  - <span class="parameter">call_doublets_thr</span> `Float`, Default: 0.25<br>
        If use_thr (previous parameter) is set to True, the threshold specified here will be used to define doublets.
    

### RNA modality Quality Control
In the following, we specify parameters used for running QC on the RNA modality data.
In the ingestion workflow we compute cell and genes QC metrics (such as % of mitochondrial genes, number of genes expressed in a cell etc.) but we do not apply filtering, as this is part of the subsequent workflow (preprocess).
Feel free to leave options blank to run with default parameters.

#### Providing a gene list
To calculate RNA QC metrics based on custom genes annotations, we need to use a gene list providing additional information on the genes expressed in the data.
Additionally, we can specify what actions we want to apply to the genes, such as what metrics to calculate.

Please visit our documentation section on [creating and using custom genes lists](../usage/gene_list_format.md) to perform quality control and visualization. 
<span class="parameter">custom_genes_file</span>`String`, Mandatory parameter, Default: resources/qc_genelist_1.0.csv<br>
    Path to the file containing the entire human gene list. Panpipes provides such a file with standard genes, and the path to this file is set as default.

##### Working with different species than human
*If working with a different species, the user must provide the appropriate gene list. For example, we offer a precompiled version of the qc gene list for mouse, the user can supply the list by specifying the path to the file as shown here:*

 `custom_genes_file:  qc_gene_list_mouse.csv`

*Find the mouse gene list in our [resources](https://github.com/DendrouLab/panpipes/blob/mouse_gene_list_upload/panpipes/resources/qc_gene_list_mouse.csv)*


It's convenient to rely on known gene lists, as this simplifies various downstream tasks, such as evaluating the percentage of mitochondrial genes in the data, identify ribosomal genes, or excluding IGG genes from HVG selection.
For the ingestion workflow, we retrieved the cell cycle genes used in `scanpy.score_genes_cell_cycle` [Satija et al. (2015), Nature Biotechnology](https://www.nature.com/articles/nbt.3192) 

| mod | feature | group  |
|-----|---------|--------|
| RNA | gene_1  | mt     |
| RNA | gene_2  | rp     |
| RNA | gene_3  | exclude|
| RNA | gene_3  | markerX|  

#### Defining actions on the genes
Next, we define "actions" on the genes as follows:

In the group column, specify what actions you want to apply to that specific gene.
For instance: `calc_proportion: mt` will calculate proportion of reads mapping to the genes whose group is "mt" in the custom genes file.

(for pipeline_ingest.py)
calc_proportions: calculate proportion of reads mapping to X genes over total number of reads, per cell
score_genes: using scanpy.score_genes function, 

(for pipeline_preprocess.py)
exclude: exclude these genes from the HVG selection, if they are deemed HV.


<span class="parameter">calc_proportions</span> `String` (comma-separated), Default: hb,mt,rp<br>
    Specify what gene proportions you want to calculate for each cell (e.g. mt for mitochondrial).
    
<span class="parameter">score_genes</span> `String`, Default: (blank) <br>
     Specify what genes should be scored.

Furthermore, there is the possibility to define a cell cycle action:

<span class="parameter">ccgenes</span> `String`, Default: default<br>
    `ccgenes` will plot the proportions of cell cycle genes, and, for each cell, determine in which cell cycle stage the respective cell is in.
    Internally, `ccgenes` uses `scanpy.tl.score_genes_cell_cycle`, which requires a file comprising cell cicle genes to be provided.
    Specify if you want to leave the default [cell cycle genes file provided by panpipes](panpipes/resources/cell_cicle_genes.tsv) (by setting this parameter to `default`) or if you want to provide your own list, in that case specify the path to that file in this parameter.
    We recommend leaving this parameter as `default`.
    If left blank, the cellcycle score will not be calculated.

### Plotting utilities for QC plots

<span class="parameter">plotqc_grouping_var</span> `String`, Default: orig.ident<br>
    Specify column in the MuData observations (MuData.obs) that stores sample information (e.g. "sample_id" or "orig.ident"). Those values will be the basis of the QC plots.
    It is also possible to use several obs columns by providing a comma separated String of multiple categorical observations (e.g. plotqc_grouping_var: sample_id,rna:channel,prot:sample)
    If left blank, the base of the plot will be the `sample_id` of your submission file.

### Plotting RNA QC metrics
All parameter values in this section should be provided as a comma separated String e.g. a,b,c.

<span class="parameter">plotqc_rna_metrics</span> `String` (comma-separated), Default: doublet_scores,pct_counts_mt,pct_counts_rp,pct_counts_hb,pct_counts_ig<br>
    What cell observations to plot for the RNA modality.
    Must be cell observations stored in .obs of the RNA AnnData.

### Plotting Protein QC metrics
Plotting QC metrices for protein data requires prot_path to be included in the submission file.
All parameter values in this section should be provided as a comma separated String e.g. a,b,c.

<span class="parameter">plotqc_prot_metrics</span> `String` (comma-separated), Default: total_counts,log1p_total_counts,n_prot_by_counts,pct_counts_isotype<br>
    By leaving this parameter as the default value, the following metrics are calculated for prot data:
    total_counts,log1p_total_counts,n_prot_by_counts,log1p_n_prot_by_counts
    If isotypes can be detected, then the following are calculated also:
    total_counts_isotype,pct_counts_isotype
    You can choose which ones you want to plot here by specifying the respective metrics.

<span class="parameter">plot_metrics_per_prot</span> `String` (comma-separated), Default: total_counts,log1p_total_counts,n_cells_by_counts,mean_counts<br>
    Since the protein antibody panels usually counts fewer features than the RNA, it might be interesting to
    visualize breakdowns of single proteins plots describing their count distribution besides other QC options.
    Specify the QC metrics you wish to plot here, such as any of the following:
    n_cells_by_counts,mean_counts,log1p_mean_counts,pct_dropout_by_counts,total_counts,log1p_total_counts

<span class="parameter">identify_isotype_outliers</span> `Boolean`, Default: True<br>
    Isotype outliers: one way to determine which cells are very sticky is to work out which cells have the most isotype UMIs associated to them.
    To label a cell as an isotype outlier, it must be in the above x% quantile by UMI counts, for at least n isotypes 
    (e.g. above 90% quantile UMIs in at least 2 isotypes).
    The actual values are specified by the following two parameters.

<span class="parameter">isotype_upper_quantile</span> `Integer`, Default: 90<br>
    See explanation for `identify_isotype_outliers`.

<span class="parameter">isotype_n_pass</span> `Integer`, Default: 2<br>
    See explanation for `identify_isotype_outliers`.

### Plot ATAC QC metrics 
We require initializing one csv file per aggregated ATAC/multiome experiment.
If you need to analyse multiple samples in the same project, aggregate them with the cellranger arc pipeline.
For multiome samples, we recommend specifying the 10X h5 input "10x_h5".
`per_barcode_metric` is only avail on cellranger arc (multiome).

<span class="parameter">is_paired</span> `Boolean`, Default: True<br>
    Are you working with only ATAC data, set to False.
    If you have multiome samples, set to True.

<span class="parameter">partner_rna</span><br>
    In case this is NOT a multiome experiment, but you have an RNA anndata that you would like to use for TSS enrichment. 
    Leave empty if no rna provided.

<span class="parameter">features_tss</span><br>
    In case this is a standalone ATAC (`is_paired`: False), please provide a feature file to run TSS enrichment. 
    Supported annotations for protein coding genes provided.

<span class="parameter">plotqc_atac_metrics</span> `String` (comma-separated), Default: n_genes_by_counts,total_counts,pct_fragments_in_peaks,atac_peak_region_fragments,atac_mitochondrial_reads,atac_TSS_fragments<br>
    Specify the ATAC metrics you want to plot and save in the metadata.
    The metrics should be provided as a comma separated string e.g. a,b,c.

### Plot Repertoire QC metrics
Repertoire data will be stored in one modality called "rep".
If you provide both TCR and BCR data, then this will be merged, nevertheless, various functions will be run on TCR and BCR separately.
Review [scirpy documentation](https://scverse.org/scirpy/latest/index.html) for specifics on the storage of the data.

<span class="parameter">ir_dist</span><br>
    Compute sequence distance metric (required for clonotype definition)
    More information on the following args are provided [here](https://scverse.org/scirpy/latest/generated/scirpy.pp.ir_dist.html#scirpy.pp.ir_dist).
    Leave blank to run with default arguments.

  - <span class="parameter">metric</span>

  - <span class="parameter">sequence</span>

<span class="parameter">clonotype_definition</span>
    Clonotype definition.
    More information on the following args are provided [here](https://scverse.org/scirpy/latest/generated/scirpy.tl.define_clonotypes.html#scirpy.tl.define_clonotypes).
    Leave blank to run with default arguments.

  - <span class="parameter">receptor_arms</span>
  - <span class="parameter">dual_ir</span>
  - <span class="parameter">within_group</span>

<span class="parameter">plotqc_rep_metrics</span><br>
Default:
- is_cell
- extra_chains
- clonal_expansion
- rep:receptor_type
- rep:receptor_subtype
- rep:chain_pairing
- rep:multi_chain

Specify which Repertoire QC metrics to plot.
Available metrics are:
- rep:clone_id_size
- rep:clonal_expansion
- rep:receptor_type
- rep:receptor_subtype
- rep:chain_pairing
- rep:multi_chain
- rep:high_confidence
- rep:is_cell
- rep:extra_chains


## Profiling Protein Ambient background
It is useful to characterize the background of your gene expression assay and antibody binding assay.
Inspect the plots and decide whether corrections such as Cellbender or SoupX (currently not included in panpipes) should be applied. 

>NOTE: This analysis can ONLY BE RUN IF YOU ARE PROVIDING RAW input starting from cellranger output, so that the "empty" droplets can be used to estimate the background.
Setting asses_background to True when you don't have RAW inputs will stop the pipeline with an error.

<span class="parameter">assess_background</span> `Boolean`, Default: False<br>
  Setting `assess_background` to True will:
  1. Create a MuData object (h5mu) from the raw data input (expected as cellranger h5 or mtx folder, if you do not have this then set to False)
  2. Plot comparative QC plots to compare distribution of UMI and feature counts in background and foreground
  3. Create heatmaps depicting the top features in the background, so that you can compare the background contamination per channel

<span class="parameter">downsample_background</span> `Boolean`, Default: True<br>
    Whether to subsample background data.
    Typically, there are many cells in the full raw cellranger output.
    Since we just want to get a picture of the background, we can subsample the data to a more reasonable size to speed up computation.
    If you want to keep all the raw data then set this parameter to False.


### Files required for profiling ambient background or running dsb normalisation
To profile ambient background or run dsb normalization, the raw_feature_bc_matrix folder from cellranger or equivalent is required.
The pipeline will automatically search for this as a .h5 or matrix folder, if the {mod}_filetype is set to "cellranger", "cellranger_multi" or "10X_h5" based on the path specified in the submission file.

If you are using a different format of file input such as a csv matrix, make sure the two files are named using the following convention:
{file_prefix}_filtered.csv" and {file_prefix}_raw.csv


### Investigate per-channel antibody staining
This can help to determine any inconsistencies in staining per channel and other QC concerns.

<span class="parameter">channel_col</span> `String`, Default: sample_id<br>
    If you want to run clr normalisation on a per-channel basis, then you need to specify which column in your submission file corresponds to the channel at the QC stage.
    It can be useful to look at the normalized data on a per-sample or per-channel basis, i.e. the 10X channel. 
    Usually, this is the sample_id column (otherwise leave the next parameter blank)

<span class="parameter">save_norm_prot_mtx</span> `Boolean`, Default: False<br>
    Set to True if you want to save the per-channel normalized values.
    It is important to note that in ingestion workflow, the per-channel normalized prot scores are not saved in the `MuData` object.
    This is because if you perform feature normalization (clr normalization margin 0 or dsb normalization) on subsets of cells, then the normalized values cannot simply be concatenated. 
    The PROT normalization is rerun on the complete object in the preprocess pipeline (or you can run this pipeline with channel_col set as None).
    Note that if you choose to run the clr on a per-channel basis, then it is not stored in the `MuData` file.


## Protein normalization

<span class="parameter">normalisation_methods</span> `String`, Default: clr, Options: dsb,clr<br>
    Choose a normalization method.
    Setting normalization method to dsb without providing raw files will stop the pipeline with an error.
    More details on this can be found [here](https://muon.readthedocs.io/en/latest/omics/citeseq.html), and more specific information on 
    [dsb here](https://muon.readthedocs.io/en/latest/api/generated/muon.prot.pp.dsb.html) and on
    [clr here](https://muon.readthedocs.io/en/latest/api/generated/muon.prot.pp.clr.html).


### Centered log ratio (CLR) normalization options

<span class="parameter">clr_margin</span> `Integer`, Default: 1<br>
    Margin determines whether to normalize per cell (as you would do for RNA normalization), 
    or by feature (recommended, due to the variable nature of prot assays). 
    CLR margin 1 is recommended for informative QC plots in this pipeline.
  - 0 = normalise row-wise (per cell)
  - 1 = normalise column-wise (per feature, recommended)

### Denoised and Scaled by Background (DSB) normalization options
 In order to run DSB you must have access to the complete raw counts, including the empty droplets from both rna and protein assays.
 See details for how to make sure your files are compatible in the _assess background_ section above.

<span class="parameter">quantile_clipping</span> `Boolean`, Default: True<br>
    Specify whether to perform quantile clipping.
    Even with normalization, some cells will have outlier values, which can be clipped as [discussed here](https://github.com/niaid/dsb).
    The maximum value will be set at the value of the 99.5% quantile, applied per feature.
    Note that this feature is in the default muon `mu.pp.dsb` code, but manually implemented here.

  
