<style>
  .parameter {
    border-top: 4px solid lightblue;
    background-color: rgba(173, 216, 230, 0.2);
    padding: 4px;
    display: inline-block;
    font-weight: bold;
  } 
</style>

# Clustering YAML  

In this documentation, the parameters of the `clustering` configuration yaml file are explained.
This file is generated running `panpipes clustering config`. <br>
The individual steps run by the pipeline are described in [clustering worlfow](https://panpipes-pipelines.readthedocs.io/en/latest/workflows/clustering.html)

When running the clustering workflow, panpipes provides a basic `pipeline.yml` file.
To run the workflow on your own data, you need to specify the parameters described below in the `pipeline.yml` file to meet the requirements of your data.

However, we do provide pre-filled versions of the `pipeline.yml` file for individual [tutorials](https://panpipes-pipelines.readthedocs.io/en/latest/tutorials/index.html).

For more information on functionalities implemented in `panpipes` to read the configuration files, such as reading blocks of parameters and reusing blocks with  `&anchors` and `*scalars`, please check [our documentation](./useful_info_on_yml.md)

You can download the different clustering pipeline.yml files here:
- Basic `pipeline.yml` file (not prefilled) that is generated when calling `panpipes clustering config`: [Download here](https://github.com/DendrouLab/panpipes/blob/main/panpipes/panpipes/pipeline_clustering/pipeline.yml)
- `pipeline.yml` for [Clustering Tutorial](https://panpipes-tutorials.readthedocs.io/en/latest/_downloads/3895aa0ba60017b15ee1aa6531dc8c25/pipeline.yml)

## Compute resources options

- <span class="parameter">resources</span><br>
Computing resources to use, specifically the number of threads used for parallel jobs.
Specified by the following three parameters:
  - <span class="parameter">threads_high</span> `Integer`, Default: 2<br>
        Number of threads used for high intensity computing tasks. 
        For each thread, there must be enough memory to load all your input files at once and create the MuData object.

  - <span class="parameter">threads_medium</span> `Integer`, Default: 2<br>
        Number of threads used for medium intensity computing tasks.
        For each thread, there must be enough memory to load your mudata and do computationally light tasks.

  - <span class="parameter">threads_low</span> `Integer`, Default: 2<br>
  	    Number of threads used for low intensity computing tasks.
        For each thread, there must be enough memory to load text files and do plotting, requires much less memory than the other two.
  - <span class="parameter">fewer_jobs</span> `Boolean`, Default: True<br>
  
  - <span class="parameter">condaenv</span> `String` (Path)<br>
    Path to conda environment that should be used to run panpipes.
    Leave blank if running native or your cluster automatically inherits the login node environment

## Loading data 
### Data format

- <span class="parameter">sample_prefix</span> `String`, Mandatory parameter, Default: mdata<br>
Prefix for the sample that comes out of the filtering/ preprocessing steps of the workflow.


- <span class="parameter">scaled_obj</span> `String`, Mandatory parameter, Default: mdata_scaled.h5mu<br>
 Path to the output file from preprocessing (e.g. `../preprocess/mdata_scaled.h5mu`).
 Ensure that the path to the file is correct.  

- <span class="parameter">full_obj</span> `String`, Default: <br>
  Speficy the full object if your scaled_obj contains only HVG.  If your scaled_obj contains all the genes then leave full_obj blank. 
  panpipes will use the full object to do marker genes analysis (rank_gene_groups) and for plotting those genes. 
- <span class="parameter">modalities</span><br>
  - <span class="parameter">rna</span> `Boolean`, Default: True<br>
  - <span class="parameter">prot</span> `Boolean`, Default: True<br>
  - <span class="parameter">atac</span> `Boolean`, Default: False<br>
  - <span class="parameter">spatial</span> `Boolean`, Default: False<br>
  Run clustering on each individual modality.

- <span class="parameter">moltimodal</span><br>
  - <span class="parameter">rna_clustering</span> `Boolean`, Default: True<br>
  - <span class="parameter">integration_method</span> `String`, Default: WNN<br>
  Options here include WNN, moda, and totalVI, and it tells us where to look for.

## Parameters for finding neighbours 

- <span class="parameter">neighors:</span> 
 Sets the number of neighbors to use when calculating the graph for clustering and umap.
  - <span class="parameter">rna:</span> 

     - <span class="parameter">use_existing </span> `Boolean`, Default: True<br>
     - <span class="parameter">dim_red </span> `String`, Default: X_pca<br>
       Defines which representation in .obsm to use for nearest neighbors
     - <span class="parameter">n_dim_red</span> `Integer`, Default: 30<br>
       Number of components to use for clustering
     - <span class="parameter">k</span> `Integer`, Default: 30<br>
       Number of neighbours
     - <span class="parameter">metric</span> `String`, Default: euclidean<br>
       Options here include euclidean and cosine
     - <span class="parameter">method</span> `String`, Default: scanpy<br>
       Options include scanpy and hnsw (from scvelo)
      
     
  - <span class="parameter">prot:</span> 

     - <span class="parameter">use_existing </span> `Boolean`, Default: True<br>
     - <span class="parameter">dim_red </span> `String`, Default: X_pca<br>
       Defines which representation in .obsm to use for nearest neighbors
     - <span class="parameter">n_dim_red</span> `Integer`, Default: 30<br>
       Number of components to use for clustering
     - <span class="parameter">k</span> `Integer`, Default: 30<br>
       Number of neighbours
     - <span class="parameter">metric</span> `String`, Default: euclidean<br>
       Options here include euclidean and cosine
     - <span class="parameter">method</span> `String`, Default: scanpy<br>
       Options include scanpy and hnsw (from scvelo)


  - <span class="parameter">atac:</span> 

     - <span class="parameter">use_existing </span> `Boolean`, Default: True<br>
     - <span class="parameter">dim_red </span> `String`, Default: X_lsi<br>
       Defines which representation in .obsm to use for nearest neighbors
     - <span class="parameter">n_dim_red</span> `Integer`, Default: 1<br>
       Number of components to use for clustering
     - <span class="parameter">k</span> `Integer`, Default: 30<br>
       Number of neighbours
     - <span class="parameter">metric</span> `String`, Default: euclidean<br>
       Options here include euclidean and cosine
     - <span class="parameter">method</span> `String`, Default: scanpy<br>
       Options include scanpy and hnsw (from scvelo)
  


  - <span class="parameter">spatial:</span> 

     - <span class="parameter">use_existing </span> `Boolean`, Default: False<br>
     - <span class="parameter">dim_red </span> `String`, Default: X_pca<br>
       Defines which representation in .obsm to use for nearest neighbors
     - <span class="parameter">n_dim_red</span> `Integer`, Default: 30<br>
       Number of components to use for clustering
     - <span class="parameter">k</span> `Integer`, Default: 30<br>
       Number of neighbours
     - <span class="parameter">metric</span> `String`, Default: euclidean<br>
       Options here include euclidean and cosine
     - <span class="parameter">method</span> `String`, Default: scanpy<br>
       Options include scanpy and hnsw (from scvelo)
  
## Parameters for umap calculation 


  - <span class="parameter">umap:</span> 

     - <span class="parameter">run </span> `Boolean`, Default: True<br>
     - <span class="parameter">rna:</span>
         - <span class="parameter">mindist </span> `Float`, Default: 0.5<br>
           Can specify an array: 0.25,0.5
      - <span class="parameter">prot:</span>
         - <span class="parameter">mindist </span> `Float`, Default: 0.5<br>
           Can specify an array: 0.25,0.5,0.8
      - <span class="parameter">atac:</span>
         - <span class="parameter">mindist </span> `Float`, Default: 0.5<br>
           Can specify an array: 0.25,0.5,0.8
      - <span class="parameter">multimodal:</span>
         - <span class="parameter">mindist </span> `Float`, Default: 0.5<br>
           Can specify an array: 0.25,0.5,0.8
      - <span class="parameter">rna:</span>
         - <span class="parameter">mindist </span> `Float`, Default: 0.5<br>
            Can specify an array: 0.25,0.5,0.8

## Parameters for clustering 

  - <span class="parameter">clusterspecs:</span>
      - <span class="parameter">rna:</span>
          - <span class="parameter">resolutions </span> `Float`, Default: 0.2, 0.6, 1<br>
           Can specify an array: 0.2,0.6,1
          - <span class="parameter">algorithm</span> `String`, Default: leiden<br>
            Options include louvain or leiden. 
      - <span class="parameter">prot:</span>
          - <span class="parameter">resolutions </span> `Float`, Default: 0.2, 0.6, 1<br>
           Can specify an array: 0.2,0.6,1
          - <span class="parameter">algorithm</span> `String`, Default: leiden<br>
            Options include louvain or leiden.

      - <span class="parameter">atac:</span>
          - <span class="parameter">resolutions </span> `Float`, Default: 0.2, 0.6, 1<br>
           Can specify an array to compute in parallel: 0.2,0.6,1
          - <span class="parameter">algorithm</span> `String`, Default: leiden<br>
            Options include louvain or leiden. 
      - <span class="parameter">multimmodal:</span>
          - <span class="parameter">resolutions </span> `Float`, Default: 0.5, 0.7<br>
           Can specify an array to compute in parallel: 0.2,0.6,1 
          - <span class="parameter">algorithm</span> `String`, Default: leiden<br>
            Options include louvain or leiden.

      - <span class="parameter">spatial:</span>
          - <span class="parameter">resolutions </span> `Float`, Default: 0.2, 0.6, 1<br>
           Can specify an array to compute in parallel: 0.2,0.6,1 
          - <span class="parameter">algorithm</span> `String`, Default: leiden<br>
            Options include louvain or leiden. 

## Parameters for finding marker genes 

In this part of the analysis we define parameters to run marker analysis. 
By default, pseudo_seurat is set to False, and we run [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html). 
When pseudo_seurat is set to True then a [python implementation](https://github.com/DendrouLab/panpipes/blob/main/panpipes/python_scripts/run_find_markers_multi.py) of `Seurat:::FindMarkers` is run

  - <span class="parameter">markerspecs:</span> <br>
     - <span class="parameter">rna:</span><br>
       - <span class="parameter">run </span> `Boolean`, Default: True<br>
       - <span class="parameter">layer </span> `String`, Default: logged_counts<br>
         Which layer stores counts for differential expression test.
       - <span class="parameter">method </span> `String`, Default: t-test_overestim_var<br>
       Options include: ‘logreg’, ‘t-test’, ‘wilcoxon’, ‘t-test_overestim_var’
       - <span class="parameter">mincels </span> `Integer`, Default: 10<br>
       If the number of clusters contains less than the number of cells maker analysis is not necessary.
       - <span class="parameter">pseudo_seurat </span> `Boolean`, Default: False<br>
       - <span class="parameter">minpct </span> `Float`, Default: 0.1<br>
       This parameter is mandatory if pseudo_seurat is set to True 
       - <span class="parameter">threshuse </span> `Float`, Default: 0.25<br>
       TThis parameter is mandatory if pseudo_seurat is set to True 
 - <span class="parameter">prot:</span><br>
   - <span class="parameter">run </span> `Boolean`, Default: True<br>
   - <span class="parameter">layer </span> `String`, Default: clr<br>
       Which layer stores counts for differential expression test.
   - <span class="parameter">mincels </span> `Integer`, Default: 10<br>
       If the number of clusters contains less than the number of cells maker analysis is not necessary.
   - <span class="parameter">method </span> `String`, Default: wilcoxon<br>
   - <span class="parameter">pseudo_seurat </span> `Boolean`, Default: False<br>
   - <span class="parameter">minpct </span> `Float`, Default: 0.1<br>
       This parameter is mandatory if pseudo_seurat is set to True 
   - <span class="parameter">threshuse </span> `Float`, Default: 0.25<br>
       This parameter is mandatory if pseudo_seurat is set to True 

 - <span class="parameter">atac:</span><br>
    - <span class="parameter">run </span> `Boolean`, Default: False<br>
    - <span class="parameter">layer </span> `String`, Default: logged_counts<br>
      Which layer stores counts for differential expression test. 
       Options include logged_counts, signac_norm , and logTF_norm,logIDF_norm
    - <span class="parameter">mincels </span> `Integer`, Default: 10<br>
       If the number of clusters contains less than the number of cells maker analysis is not necessary.
    - <span class="parameter">method </span> `String`, Default: wilcoxon<br>
        Options include: ‘logreg’, ‘t-test’, ‘wilcoxon’, ‘t-test_overestim_var’
    - <span class="parameter">pseudo_seurat </span> `Boolean`, Default: False<br>
    - <span class="parameter">minpct </span> `Float`, Default: 0.1<br>
       This parameter is mandatory if pseudo_seurat is set to True 
    - <span class="parameter">threshuse </span> `Float`, Default: 0.25<br>
       This parameter is mandatory if pseudo_seurat is set to True 


 - <span class="parameter">multimodal:</span><br>
   - <span class="parameter">mincels </span> `Integer`, Default:10<br>
       If the number of clusters contains less than the number of cells maker analysis is not necessary.
    - <span class="parameter">method </span> `String`, Default: wilcoxon<br>
        Options include: ‘logreg’, ‘t-test’, ‘wilcoxon’, ‘t-test_overestim_var’
    - <span class="parameter">pseudo_seurat </span> `Boolean`, Default: False<br>
    - <span class="parameter">minpct </span> `Float`, Default: 0.1<br>
       This parameter is mandatory if pseudo_seurat is set to True 
    - <span class="parameter">threshuse </span> `Float`, Default: 0.25<br>
       This parameter is mandatory if pseudo_seurat is set to True


 - <span class="parameter">spatial:</span><br>
   - <span class="parameter">run </span> `Boolean`, Default: True<br>
   - <span class="parameter">layer </span> `String`, Default: norm_pearson_resid<br>
       Options include logged_counts, signac_norm , and logTF_norm,logIDF_norm
   - <span class="parameter">method </span> `String`, Default: t-test_overestim_var<br>
        Options include: ‘logreg’, ‘t-test’, ‘wilcoxon’, ‘t-test_overestim_var’
   - <span class="parameter">mincels </span> `Integer`, Default: 10<br>
       If the number of clusters contains less than the number of cells maker analysis is not necessary.
   - <span class="parameter">pseudo_seurat </span> `Boolean`, Default: False<br>
   - <span class="parameter">minpct </span> `Float`, Default: 0.1<br>
      This parameter is mandatory if pseudo_seurat is set to True 
   - <span class="parameter">threshuse </span> `Float`, Default: 0.25<br>
       This parameter is mandatory if pseudo_seurat is set to True 
## Plot specifications
Used to define which metadata columns are used in the visualisations 
 - <span class="parameter">plotspecs:</span><br>
   - <span class="parameter">layers: </span><br>
     - <span class="parameter">rna </span> `String`, Default: logged_counts<br>
     - <span class="parameter">prot </span> `String`, Default: clr<br>
     - <span class="parameter">atac </span> `String`, Default: signac_norm<br>
     - <span class="parameter">spatial </span> `?`, Default: ?<br> 
  - <span class="parameter">top_n_markers </span> `Integer`, Default: 10<br>
