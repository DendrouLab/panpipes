Using custom genes annotations
------------------------------

It's often practical to rely on known gene lists, for a series of tasks, like evaluating % of mitochondrial genes or
ribosomal genes, or excluding genes from HVG selection such as those constituting the IG chains. 
We collect some of these genes in example input files:

- the cellcycle genes used in scanpy.score_genes_cell_cycle [Satija et al. (2015), Nature Biotechnology.] 
are stored in  panpipes/resources/cell_cicle_genes.tsv

| gene_name | cc_phase |
| --------- | -------- |
| MCM5      | s        |
| PCNA      | s        |
| TYMS      | s        |
| CCNB2     | g2m      |
| CKAP2L    | g2m      |
| ...       | ...      |


- an example of a single file storing useful gene lists in panpipes/resources/qc_genelist_1.0.csv 


| mod | feature | group   |
| --- | ------- | ------- |
| RNA | gene_1  | mt      |
| RNA | gene_2  | rp      |
| RNA | gene_1  | exclude |
| RNA | gene_1  | markerX |
| ... | ...     | ...     |



## Explaining actions

Panpipes uses "actions" to define which tasks to use which gene list for.
Specify the "group" name of the genes you want to use to apply the action i.e. calc_proportion: mt will calculate
proportion of reads mapping to the genes whose group is "mt"


(for pipeline_ingest.py)
**calc_proportions:** calculate proportion of reads mapping to X genes over total number of reads, per cell
**score_genes:** using scanpy.tl.score_genes function, the average expression of a set of genes, subtracted of the average expression of a reference set of genes. First introduced in Satija et al. Nature Biotechnology (2015).

(for pipeline_integration.py)
**exclude:** exclude these genes from the HVG selection, if they are deemed Highly Variable.
**plot_markers:** plot these genes

If left blank, these actions will not be performed (i.e. no calculation of % of mt genes per cell will be included in the ingestion of the data)

### Cell cycle action

**ccgenes:**  
Setting the `ccgenes` to `default` will calculate the phase of the cell cycle in which the cell is by using `scanpy.tl.score_genes_cell_cycle` using the file provided in panpipes/resources/cell_cicle_genes.tsv

Users can create their own list and specify the path in the `ccgenes` param to score the cells with a custom list.
If left blank, the cellcycle score will not be calculated.



