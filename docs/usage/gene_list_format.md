Using custom genes annotations: gene list formats
=================

It's often practical to rely on known gene lists, for a series of tasks, like evaluating % of mitochondrial genes or
ribosomal genes, or excluding genes from HVG selection such as those constituting the IG chains. 

### Custom gene lists

We provide an example of a preformatted gene lists file in [resources/qc_genelist_1.0.csv](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/qc_genelist_1.0.csv).

All <sup>[1](#footnote1)</sup> files provided to the pipeline should be in a 3 columns format, where the column headers are "mod" (modality: "rna", "prot", or "atac"), feature and group. The group column is used to distinguish different gene groups.


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


The gene lists are not pre-determined within panpipes in order to maximise
flexibility as all organisms will require separate lists, but there are example lists provided on our [github page](https://github.com/DendrouLab/panpipes/tree/main/panpipes/resources)

### Cell cycle genes

The cellcycle genes used in [scanpy.score_genes_cell_cycle](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes_cell_cycle.html) 
are stored in [resources/cell_cycle_genes.csv](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/cell_cycle_genes.tsv)

Differently from the other custom gene file, the cell cycle file should be a **tab separated file with two columns**:

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

## Explaining actions

Panpipes uses "actions" to define which tasks to use which gene list for.
Specify the "group" name of the genes you want to use to apply the action i.e. calc_proportion: mt will calculate
proportion of reads mapping to the genes whose group is "mt"


The genes are scored for each modality using 

- [scanpy.pp.calculate_qc_metrics](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics).
For example, for the rna modality, including a list of mitochondiral
genes in the group `mt`, will add `pct_counts_mt`, and `total_counts_mt`
to the `mdata["rna"].obs` assay.

- [scanpy.tl.score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html).


(for pipeline_ingest.py)
**calc_proportions:** calculate proportion of reads mapping to X genes over total number of reads, per cell
**score_genes:** using scanpy.tl.score_genes function, the average expression of a set of genes, subtracted of the average expression of a reference set of genes. First introduced in Satija et al. Nature Biotechnology (2015).

(for pipeline_preprocess.py)
**exclude:** exclude these genes from the HVG selection, if they are deemed Highly Variable.

For the exclude action, if set to `default` the workflow will look for genes whose group is set to `exclude` in the supplied qc_genelist file. Alternatively, if you are specifying your custom gene list and you want to exclude another set of genes, for example a group you call `TCR_genes`, specify this group (i.e. `exclude: TCR_genes`)

If left blank, these actions will not be performed (i.e. no calculation of % of mt genes per cell will be included in the ingestion of the data)

### Cell cycle action

**ccgenes:**  
Setting the `ccgenes` param to `default` in the ingest workflow will calculate the phase of the cell cycle in which the cell is by using `scanpy.tl.score_genes_cell_cycle` using the file provided in panpipes/resources/cell_cicle_genes.tsv

Users can create their own list, and need to specify the path to this new file in in the `ccgenes` param to score the cells with their custom list.

If left blank, the cellcycle score will not be calculated.


## Supplying custom gene lists

The custom genelist file can be supplied by the user in three worflows: 

Ingest workflow
-------------

pipeline_ingest config file: (pipeline.yml)

```
custom_genes_file: resources/qc_genelist_1.0.csv
```

Preprocess workflow
----------------

pipeline_preprocess config file: (pipeline.yml)

```
exclude_file: resources/qc_genelist_1.0.csv
```

*Note that we have formatted an example file containing all genes to use in both workflows, and therefore supply the same file to both workflows but users can have independent files for each of them.*

Vizualization workflow
---------------

pipeline_vis config file: (pipeline.yml)

```
# the full list will be plotted in dot plots and matrix plots, one plot per group
full:
 - long_file1.csv
 - long_file2.csv
# the shorter list will be plotted on umaps as well as other plot types, one plot per group
minimal:
 - short_file1.csv

```

These require the same 3-column format as above.
Generally in the visualisation pipeline all gene groups in the input are plotted. In heatmaps and dot
plots, one dotplot per group is plotted. For umaps, one plot per gene is
plotted, and a new file is saved per group.

## Final notes

Be deliberate and informative with the choice of group names for any gene set use, since the .obs column generated as output will be named based on the group of the gene list input file. 
The columns added to the mudata object will be used for filtering in the `preprocess` workflow (see [filtering instructions](./filter_dict_instructions.md) ).

If the mitochondrial genes are in group "mt" as in the example given in the resource file, then the column generated with the **calc_proportions** action and containing the percentage of MT genes will be named "pct_counts_mt".

So specifying *pct_counts_mito* instead of *pct_counts_mt* will not filter the mudata based on mitochondrial %, because the workflow can't find the supplied column in the data.

