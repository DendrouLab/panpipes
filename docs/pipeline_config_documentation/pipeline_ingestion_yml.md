<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  }
</style>

# pipeline.yml for ingestion workflow

## Compute resources options

* <p class="parameter">resources</p>
  
    Computing resources to use, specifically the number of threads used for parallel jobs.
    Specified by the following three parameters:

  - <p class="parameter">threads_high</p> Integer, Default: 1
        <p>Number of threads used for high intensity computing tasks. 
        For each thread, there must be enough memory to load all your input files at once and create the MuData object.
        </p><br>
  
  - <p class="parameter">threads_medium</p> Integer, Default: 1
        <p>Number of threads used for medium intensity computing tasks.
        For each thread, there must be enough memory to load your mudata and do computationally light tasks.
        </p><br>
  
  - <p class="parameter">threads_low</p> Integer, Default: 1
  	    <p>Number of threads used for low intensity computing tasks.
        For each thread, there must be enough memory to load text files and do plotting, requires much less memory than the other two.
        </p><br>

* <p class="parameter">condaenv</p> String
  
    Path to conda environment that should be used to run panpipes.
    Leave blank if running native or your cluster automatically inherits the login node environment


## Loading and merging data options
### Project name and data format

* <p class="parameter">project</p> String, Default: "test"
    
    Project name.


* <p class="parameter">sample_prefix</p> String, Default: "test"
    
    Prefix for sample names.


* <p class="parameter">use_existing_h5mu</p> Boolean (True/False), Default: False
    
    If you want to read an existing MuData object (.h5mu file), set to True. Move the object to the folder where you intend to run the ingestion workflow, and call it ${sample_prefix}_unfilt.h5mu where ${sample_prefix} = sample_prefix argument above.


* <p class="parameter">submission_file</p> String, Mandatory parameter
    
    Submission file name (e.g. sample_file_qc.txt).
    Ensure that the submission file must be in the right format.
    For ingest, the submission file must contain the follwing columns: sample_id  rna_path  rna_filetype  (prot_path  prot_filetype tcr_path  tcr_filetype etc.)
    An example submission file can be found at resources/sample_file_mm.txt


* <p class="parameter">metadatacols</p> String (comma-separated)
    
    Specify the metadata columns from the submission file that you want to include in for the downstream processing of the data.
    The specified columns will we included in the created AnnData object that stores all data. 
    Provide the columns as a comma-separated String, for example: batch,disease,sex


* <p class="parameter">concat_join_type</p> String, Default: inner
    
    Specifies which concat join is performed on the MuData objects.
    We recommended `inner`.
    See the [AnnData documentation](https://anndata.readthedocs.io/en/latest/concatenation.html#inner-and-outer-joins) for details.


### Modalities in the project

It is crucial to specify which modalities are present in the data and need to be processed by panpipes.
The Quality Control (QC) is run for each modality independently.

* <p class="parameter">modalities:</p>
    
    Specify which modalities are included in the data by setting the respective modality to True.
    Leave empty (None) or False to signal this modality is not part of the experiment.
    The modalities are processed in the order of the following list:

    * <p class="parameter">rna</p> Boolean, Default: False
    * <p class="parameter">prot</p> Boolean, Default: False
    * <p class="parameter">bcr</p> Boolean, Default: False
    * <p class="parameter">tcr</p> Boolean, Default: False
    * <p class="parameter">atac</p> Boolean, Default: False
    
### Integrating barcode level data

* <p class="parameter">barcode_mtd</p>
  
    Optionally, you can ingest cell-barcode-level metadata, such as demultiplexing with hashtags, chemical tags or lipid tagging.
    In case you have cell level metadata such as results from a demultiplexing algorithm, you can incorporate it into the MuData object.
    The demultiplexing results should be stored in one csv file containing 2 columns, barcode_id, and sample_id.
    Note that the sample_id should match the sample_id column in the submission file.
  
    * <p class="parameter">include</p> Boolean, Default: False<br>
        <p>Set to True if you want to include barcode level data.
        </p><br>
      
    * <p class="parameter">path</p> String
        <p>Insert path to the csv file containing the demultiplexing results.
           Requires include to be set to True.
        </p><br>
      
    * <p class="parameter">metadatacols</p> String
        <p>Provide the metadata columns in the file specified via path that should be incorporated in the MuData object.
        </p><br>
    

### Loading Protein data - additional options
Optionally, you can provide additional information on the Protein modality by specifying the following parameters.

* <p class="parameter">protein_metadata_table</p> String
    
    File containing additional information on the proteins.
    By default, the ingestion workflow will use the first column of the cellranger features.tsv.gz.
    If you want to add additional information about the antibodies (e.g. whether they are hashing antibodies or isotypes), you need to provide an additional table containing these information.
    The first column of the table must be equivalent to the first column of cellrangers features.tsv.gz and it must be stored as a txt file.
    For example, you could create a table that include a column called isotype (containing True or False) to run QC isotypes, and a column called hashing_ab (containing True or False) if your dataset contains hashing antibodies which you wish to split into a separate modality.


* <p class="parameter">index_col_choice</p> String
    
    Column name from protein_metadata_table that should be included in the MuData object.
    If you want to update the MuData object storing all our data with a column from the table specified in protein_metadata_table, specify the column name here.
    Ensure that there are no overlaps with the gene symbols used for the rna modality, otherwise you might face issues down the line.
    If there are overlaps with the rna gene symbols, then you might have trouble down the line
    Consider adding e.g. a prefix to make the protein labels unique (e.g. "prot_CD4").
 

* <p class="parameter">load_prot_from_raw</p> Boolean, Default: False
    
    In case separate prot and rna 10X outputs are provided, the pipeline will load the filtered rna 10X outputs and the raw prot 10X counts.

    
* <p class="parameter">subset_prot_barcodes_to_rna</p> Default: False
    
    If providing separate prot and rna 10X output (`load_prot_from_raw`: True), consider setting `subset_prot_barcodes_to_rna` to True.
    Set to True if you want to treat filtered rna barcodes as the "real" cells (which we assume one wants to do in most cases), and, additionally, the out prot assay barcodes should match the rna ones.
    If set to False, the full prot matrix will be loaded into the MuData object.


## Quality Control (QC) options
### Processing of 10X cellranger files

* <p class="parameter">plot_10X_metrics</p> Boolean, Default: False<br>
    
    Specify if you want to plot 10x cellranger QC metrics.
    Set this parameter to False if you're not using cellranger folders as input.
    If starting from cellranger output, you can parse multiple metrics_summary.csv files and plot them together.
    The ingestion workflow will search for the metrics_summary file(s) in the "outs" folder and will stop if it doesn't find it.
 

### Doublet detection on RNA modality
For doublet detection, we use the [Scrublet tool](https://www.sciencedirect.com/science/article/pii/S2405471218304745?via%3Dihub) for detection of doublets in the data. Doublets are cells that contain two distinct cell barcodes, and are often a result of two cells being encapsulated in the same droplet during the library preparation step.

* <p class="parameter">scr</p><br>
    
    Specify the following parameters for doublet detection using the tool Scrublet.

    * <p class="parameter">expected_doublet_rate</p> Float, Default: 0.06
        
        Fraction of observations (cells) expected to be doublets.
        Note that results are not particularly sensitive to this parameter. 
  
    * <p class="parameter">sim_doublet_ratio</p> Float, Default: 2
        
        Number of doublets to simulate, relative to the number of total observations.
        Note that setting this too high drastically increases computationally complexity.
        The minimum value we tested is 0.5.
      
    * <p class="parameter">n_neighbours</p>Integer, Default: 20
        
        Number of neighbors used to construct the k-nearest-neighbor (KNN) classifier for the observations and simulated doublets.
        The value (round(0.5*sqrt(n_cells))) usually works well.
    
    * <p class="parameter">min_counts</p>Integer, Default: 2
        
        Used for gene filtering prior to PCA.
        Genes expressed less than `min_counts` times in `min_cells` cells (see next parameter) are excluded from the data for further processing.

    * <p class="parameter">min_cells</p>Integer, Default: 3
      
        Used for gene filtering prior to PCA.
        Genes expressed less than `min_counts` (see previous parameter) times in `min_cells` cells are excluded from the data for further processing."
  
    * <p class="parameter">min_gene_variability_pctl</p>Integer, Default: 85
      
        Used for gene filtering prior to PCA.
        Keeps only the highly variable genes in the top min_gene_variability_pctl percentile for downstream processing, as measured by the v-statistic [Klein et al., Cell 2015].
  
    * <p class="parameter">n_prin_comps</p> Integer, Default: 30
      
        Number of principal components used to embed the observations in a lower-dimensional space.
        PCA is run prior to k-nearest-neighbor graph construction.
  
    * <p class="parameter">use_thr</p> Boolean, Default: True
      
        Specify if you want to use a user-defined threshold in plots to distinguish between true and false doublets.
        The actual threshold is defined in the following parameter.
        If set to False, doublet detection will be run with the default threshold calculated by Scrublet.
        Note that this threshold applies to plots only, meaning no actual filtering takes place here.
  
    * <p class="parameter">call_doublets_thr</p> Default: 0.25
      
        If use_thr (previous parameter) is set to True, the threshold specified here will be used to define doublets.
    

### RNA modality Quality Control
In the following, we specify parameters used for running QC on the RNA modality data.
In the ingestion workflow we compute cell and genes QC metrics (such as % of mitochondrial genes, number of genes expressed in a cell etc.) but we do not apply filtering, as this is part of the subsequent workflow (preprocess).
Feel free to leave options blank to run with default parameters.

#### Providing a gene list
To calculate RNA QC metrics, we need to define a gene list providing additional information on the genes in the data.
Additionally, we can specify what actions we want to apply to the genes, such as what metrics to calculate.

* <p class="parameter">custom_genes_file</p>String, Default: resources/qc_genelist_1.0.csv
    
    Path to the file containing the entire gene list. Panpipes provides such a file with standard genes, and the path to this file is set as default.
 

Usually, it's convenient to rely on known gene lists, as this simplifies various downstream tasks, such as evaluating the percentage of mitochondrial genes in the data, identify ribosomal genes, or excluding IGG genes from HVG selection.
For the ingestion workflow, we retrieved the cell cycle genes used in `scanpy.score_genes_cell_cycle` [Satija et al. (2015), Nature Biotechnology](https://www.nature.com/articles/nbt.3192) and stored them in a file: panpipes/resources/cell_cicle_genes.tsv.
Additionally, we also provide an example for an entire gene list: panpipes/resources/qc_genelist_1.0.csv 

| mod | feature | group  |
|-----|---------|--------|
| RNA | gene_1  | mt     |
| RNA | gene_2  | rp     |
| RNA | gene_3  | exclude|
| RNA | gene_3  | markerX|  

#### Defining actions on the genes
Next, we define "actions" on the genes as follows:

In the group column, specify what actions you want to apply to that specific gene.
For instance: calc_proportion: mt will calculate proportion of reads mapping to the genes whose group is "mt".

(for pipeline_ingest.py)
calc_proportions: calculate proportion of reads mapping to X genes over total number of reads, per cell
score_genes: using scanpy.score_genes function, 

(for pipeline_preprocess.py)
exclude: exclude these genes from the HVG selection, if they are deemed HV.


* <p class="parameter">calc_proportions</p> Default: hb,mt,rp
    
    Specify what gene proportions you want to calculate for each cell (e.g. mt for mitochondrial).
    

* <p class="parameter">score_genes</p> Default: MarkersNeutro
    
    Specify what genes should be scored.

Furthermore, there is the possibility to define a cell cycle action:

* <p class="parameter">ccgenes</p> String, Default: default
    
    `ccgenes` will plot the proportions of cell cycle genes, and, for each cell, determine in which cell cycle stage the respective cell is in.
    Internally, `ccgenes` uses `scanpy.tl.score_genes_cell_cycle`, which requires a file comprising cell cicle genes to be provided.
    Specify if you want to leave the default [cell cycle genes file provided by panpipes](panpipes/resources/cell_cicle_genes.tsv) (by setting this parameter to `default`) or if you want to provide your own list, in that case specify the path to that file in this parameter.
    We recommend leaving this parameter as `default`.
    If left blank, the cellcycle score will not be calculated.
 

### Plotting utilities for QC plots

* <p class="parameter">plotqc_grouping_var</p> String, Default: orig.ident
    
    Specify column in the MuData observations (MuData.obs) that stores sample information (e.g. "sample_id" or "orig.ident"). Those values will be the basis of the QC plots.
    It is also possible to use several obs columns by providing a comma separated String of multiple categorical observations (e.g. plotqc_grouping_var: sample_id,rna:channel,prot:sample)
    If left blank, the base of the plot will be the `sample_id` of your submission file.
    

### Plotting RNA QC metrics
All parameter values in this section should be provided as a comma separated String e.g. a,b,c.

* <p class="parameter">plotqc_rna_metrics</p> String (comma-separated), Default: doublet_scores,pct_counts_mt,pct_counts_rp,pct_counts_hb,pct_counts_ig
    
    What cell observations to plot for the RNA modality.
    Must be cell observations stored in .obs of the RNA AnnData.
  

### Plotting Protein QC metrics
Plotting QC metrices for protein data requires prot_path to be included in the submission file.
All parameter values in this section should be provided as a comma separated String e.g. a,b,c.

* <p class="parameter">plotqc_prot_metrics</p> String (comma-separated), Default: total_counts,log1p_total_counts,n_prot_by_counts,pct_counts_isotype
    
    By leaving this parameter as the default value, the following metrics are calculated for prot data:
    total_counts,log1p_total_counts,n_prot_by_counts,log1p_n_prot_by_counts
    If isotypes can be detected, then the following are calculated also:
    total_counts_isotype,pct_counts_isotype
    You can choose which ones you want to plot here by specifying the respective metrics.
    

* <p class="parameter">plot_metrics_per_prot</p> String (comma-separated), Default: total_counts,log1p_total_counts,n_cells_by_counts,mean_counts
    
    Since the protein antibody panels usually counts fewer features than the RNA, it might be interesting to
    visualize breakdowns of single proteins plots describing their count distribution besides other QC options.
    Specify the QC metrics you wish to plot here, such as any of the following:
    n_cells_by_counts,mean_counts,log1p_mean_counts,pct_dropout_by_counts,total_counts,log1p_total_counts
     

* <p class="parameter">identify_isotype_outliers</p> Boolean, Default: True
    
    Isotype outliers: one way to determine which cells are very sticky is to work out which cells have the most isotype UMIs associated to them.
    To label a cell as an isotype outlier, it must be in the above x% quantile by UMI counts, for at least n isotypes 
    (e.g. above 90% quantile UMIs in at least 2 isotypes).
    The actual values are specified by the following two parameters.
    

* <p class="parameter">isotype_upper_quantile</p> Integer, Default: 90
    
  See explanation for `identify_isotype_outliers`.

    
* <p class="parameter">isotype_n_pass</p> Integer, Default: 2<br>
    
    See explanation for `identify_isotype_outliers`.


### Plot ATAC QC metrics 
We require initializing one csv file per aggregated ATAC/multiome experiment.
If you need to analyse multiple samples in the same project, aggregate them with the cellranger arc pipeline.
For multiome samples, we recommend specifying the 10X h5 input "10x_h5".
`per_barcode_metric` is only avail on cellranger arc (multiome).

* <p class="parameter">is_paired</p> Boolean, Default: True
  
    Are you working with only ATAC data, set to False.
    If you have multiome samples, set to True.


* <p class="parameter">partner_rna</p>
  
    In case this is NOT a multiome experiment, but you have an RNA anndata that you would like to use for TSS enrichment. 
    Leave empty if no rna provided.


* <p class="parameter">features_tss</p>
    
    In case this is a standalone ATAC (`is_paired`: False), please provide a feature file to run TSS enrichment. 
    Supported annotations for protein coding genes provided.


* <p class="parameter">plotqc_atac_metrics</p> String (comma-separated), Default: n_genes_by_counts,total_counts,pct_fragments_in_peaks,atac_peak_region_fragments,atac_mitochondrial_reads,atac_TSS_fragments
    
    Specify the ATAC metrics you want to plot and save in the metadata.
    The metrics should be provided as a comma separated string e.g. a,b,c.


### Plot Repertoire QC metrics
Repertoire data will be stored in one modality called "rep".
If you provide both TCR and BCR data, then this will be merged, nevertheless, various functions will be run on TCR and BCR separately.
Review [scirpy documentation](https://scverse.org/scirpy/latest/index.html) for specifics on the storage of the data.

* <p class="parameter">ir_dist</p>
  
    Compute sequence distance metric (required for clonotype definition)
    More information on the following args are provided [here](https://scverse.org/scirpy/latest/generated/scirpy.pp.ir_dist.html#scirpy.pp.ir_dist).
    Leave blank to run with default arguments.

    * <p class="parameter">metric</p>
  
    * <p class="parameter">sequence</p>

  
* <p class="parameter">clonotype_definition</p>
    
    Clonotype definition.
    More information on the following args are provided [here](https://scverse.org/scirpy/latest/generated/scirpy.tl.define_clonotypes.html#scirpy.tl.define_clonotypes).
    Leave blank to run with default arguments.

    * <p class="parameter">receptor_arms</p>
    
    * <p class="parameter">dual_ir</p>
    
    * <p class="parameter">within_group</p>

  
* <p class="parameter">plotqc_rep_metrics</p>
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
    * rep:clone_id_size
    * rep:clonal_expansion
    * rep:receptor_type
    * rep:receptor_subtype
    * rep:chain_pairing
    * rep:multi_chain
    * rep:high_confidence
    * rep:is_cell
    * rep:extra_chains

    
    
## Profiling Protein Ambient background
It is useful to characterize the background of your gene expression assay and antibody binding assay.
Inspect the plots and decide whether corrections such as Cellbender or SoupX (currently not included in panpipes) should be applied. 

>NOTE: This analysis can ONLY BE RUN IF YOU ARE PROVIDING RAW input starting from cellranger output, so that the "empty" droplets can be used to estimate the background.
Setting asses_background to True when you don't have RAW inputs will stop the pipeline with an error.

* <p class="parameter">assess_background</p> Boolean, Default: False

    Setting `assess_background` to True will:
  
    1. Create a MuData object (h5mu) from the raw data input (expected as cellranger h5 or mtx folder, if you do not have this then set to False)
    2. Plot comparative QC plots to compare distribution of UMI and feature counts in background and foreground
    3. Create heatmaps depicting the top features in the background, so that you can compare the background contamination per channel


* <p class="parameter">downsample_background</p> Boolean, Default: True
    
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

* <p class="parameter">channel_col</p> String, Default: sample_id
  
    If you want to run clr normalisation on a per-channel basis, then you need to specify which column in your submission file corresponds to the channel at the QC stage.
    It can be useful to look at the normalized data on a per-sample or per-channel basis, i.e. the 10X channel. 
    Usually, this is the sample_id column (otherwise leave the next parameter blank)


* <p class="parameter">save_norm_prot_mtx</p> Boolean, Default: False
  
    Set to True if you want to save the per-channel normalized values.
    It is important to note that in ingestion workflow, the per-channel normalized prot scores are not saved in the `MuData` object.
    This is because if you perform feature normalization (clr normalization margin 0 or dsb normalization) on subsets of cells, then the normalized values cannot simply be concatenated. 
    The PROT normalization is rerun on the complete object in the preprocess pipeline (or you can run this pipeline with channel_col set as None).
    Note that if you choose to run the clr on a per-channel basis, then it is not stored in the `MuData` file.


## Protein normalization

* <p class="parameter">normalisation_methods</p> String, Default: clr, Options: dsb,clr
    
    Choose a normalization method.
    Setting normalization method to dsb without providing raw files will stop the pipeline with an error.
    More details on this can be found [here](https://muon.readthedocs.io/en/latest/omics/citeseq.html), and more specific information on 
    [dsb here](https://muon.readthedocs.io/en/latest/api/generated/muon.prot.pp.dsb.html) and on
    [clr here](https://muon.readthedocs.io/en/latest/api/generated/muon.prot.pp.clr.html).


### Centered log ratio (CLR) normalization options

* <p class="parameter">clr_margin</p> Integer, Default: 0
  
    Margin determines whether to normalize per cell (as you would do for RNA normalization), 
    or by feature (recommended, due to the variable nature of prot assays). 
    CLR margin 0 is recommended for informative QC plots in this pipeline.
    - 0 = normalise rowwise (per feature, recommended)
    - 1 = normalise colwise (per cell)

### Denoised and Scaled by Background (DSB) normalization options
 In order to run DSB you must have access to the complete raw counts, including the empty droplets from both rna and protein assays.
 See details for how to make sure your files are compatible in the _assess background_ section above.

* <p class="parameter">quantile_clipping</p> Boolean, Default: True
  
    Specify whether to perform quantile clipping.
    Even with normalization, some cells will have outlier values, which can be clipped as [discussed here](https://github.com/niaid/dsb).
    The maximum value will be set at the value of the 99.5% quantile, applied per feature.
    Note that this feature is in the default muon `mu.pp.dsb` code, but manually implemented here.

  