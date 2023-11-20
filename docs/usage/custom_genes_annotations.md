Using custom genes annotations
------------------------------

It's often practical to rely on known gene lists, for a series of tasks, like evaluating % of mitochondrial genes or
ribosomal genes, or excluding genes from HVG selection such as those constituting the IG chains. 
We collect some of these genes in example input files:

- the cellcycle genes used in scanpy.score_genes_cell_cycle [Satija et al. (2015), Nature Biotechnology.] 
are stored in  panpipes/resources/cell_cicle_genes.tsv

the file should be a tab separated file with two columns:

**gene_name**:  the name of the gene
**cc_phase**: which phase of the cell cycle is the gene expression indicative of. 

| gene_name | cc_phase |
| --------- | -------- |
| MCM5      | s        |
| PCNA      | s        |
| TYMS      | s        |
| CCNB2     | g2m      |
| CKAP2L    | g2m      |
| ...       | ...      |


- an example of a single file storing useful gene lists in panpipes/resources/qc_genelist_1.0.csv 

The file should be formatted as a comma separated file, with 3 columns: 

**mod**: the modality for the feature in use 
**feature**: feature name, i.e. gene
**group**: the group the gene belongs to

| mod | feature | group   |
| --- | ------- | ------- |
| RNA | gene_1  | mt      |
| RNA | gene_2  | rp      |
| RNA | gene_1  | exclude |
| RNA | gene_1  | markerX |
| ... | ...     | ...     |


The custom genelist file should be supplied by the user in two worflows: 

pipeline_ingest config file: (pipeline.yml)

```
custom_genes_file: resources/qc_genelist_1.0.csv
```


pipeline_preprocess config file: (pipeline.yml)

```
exclude_file: resources/qc_genelist_1.0.csv
```

We have formatted an example file and therefore supply the same file to both workflows but users can have independent files for each of them.

## Explaining actions

Panpipes uses "actions" to define which tasks to use which gene list for.
Specify the "group" name of the genes you want to use to apply the action i.e. calc_proportion: mt will calculate
proportion of reads mapping to the genes whose group is "mt"


(for pipeline_ingest.py)
**calc_proportions:** calculate proportion of reads mapping to X genes over total number of reads, per cell
**score_genes:** using scanpy.tl.score_genes function, the average expression of a set of genes, subtracted of the average expression of a reference set of genes. First introduced in Satija et al. Nature Biotechnology (2015).

(for pipeline_preprocess.py)
**exclude:** exclude these genes from the HVG selection, if they are deemed Highly Variable.

For the exclude action, if set to `default` the workflow will look for genes whose group is set to `exclude` in the supplied qc_genelist file. Alternatively, if you are specifying your custom gene list and you want to exclude another set of genes, for example a group you call `TCR_genes`, specify this group (i.e. `exclude: TCR_genes`)

If left blank, these actions will not be performed (i.e. no calculation of % of mt genes per cell will be included in the ingestion of the data)

### Cell cycle action

**ccgenes:**  
Setting the `ccgenes` to `default` will calculate the phase of the cell cycle in which the cell is by using `scanpy.tl.score_genes_cell_cycle` using the file provided in panpipes/resources/cell_cicle_genes.tsv

Users can create their own list, and need to specify the path to this new file in in the `ccgenes` param to score the cells with their custom list.

If left blank, the cellcycle score will not be calculated.



