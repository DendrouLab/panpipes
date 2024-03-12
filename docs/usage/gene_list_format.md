Using custom genes annotations: gene list formats
=================

It's often practical to rely on known gene lists, for a series of tasks, like evaluating % of mitochondrial genes or
ribosomal genes, or excluding genes from HVG selection such as those constituting the IG chains. 

### Custom gene lists

We provide an example of a preformatted gene lists file in [resources/qc_genelist_1.0.csv](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/qc_genelist_1.0.csv).

All Custom Gene Lists files provided to the pipeline should be in a 3 columns format, where the column headers are "mod" (modality: "rna", "prot", or "atac"), feature and group. The group column is used to distinguish different gene groups.

- **mod**: the modality for the feature in use. Modalities are always specified in lowercase.
- **feature**: feature name, i.e. a gene or a protein id. 
- **group**: the group the gene belongs to. Group can be upper or lowercase or a mix of both and will be interpreted as a string.

| mod | feature | group   |
| --- | ------- | ------- |
| rna | gene_1  | mt      |
| rna | gene_2  | rp      |
| rna | gene_1  | exclude |
| rna | gene_1  | markerX |
| ... | ...     | ...     |

Users can provide any gene format in the custom gene files, depending on the reference they used to produce their count matrices. 
Gene ids can come in a variety of formats, including upper or lower cases, letter,numbers or a combination of both. 

```
GeneCards Symbol: RHO
HGNC: 10012 
NCBI Gene: 6010 
Ensembl: ENSG00000163914
```

Therefore, the gene lists are not pre-determined within panpipes in order to maximise flexibility and users should provide their own lists.

For a typical usecase, we provide example lists on our [github page](https://github.com/DendrouLab/panpipes/tree/main/panpipes/resources) which are also used by default as specified in the [next sections](#explaining-custom-gene-lists-actions).



### Cell cycle genes

The human-only cellcycle genes used in [scanpy.score_genes_cell_cycle](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes_cell_cycle.html) 
are stored in [resources/cell_cycle_genes.csv](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/cell_cycle_genes.tsv)

However, if the data is mouse only then the cellcycle gene list can be found in [resources/mouse_cell_cycle_genes.tsv](https://github.com/DendrouLab/panpipes/blob/mouse_cell_cycle/panpipes/resources/mouse_cell_cycle_genes.tsv)

Differently from the other custom gene file, the cell cycle file should be a **tab separated file with two columns**:

- **gene_name**:  the name of the gene
- **cc_phase**: which phase of the cell cycle is the gene expression indicative of.

| gene_name | cc_phase |
| --------- | -------- |
| MCM5      | s        |
| PCNA      | s        |
| TYMS      | s        |
| CCNB2     | g2m      |
| CKAP2L    | g2m      |
| ...       | ...      |

## Using custom gene lists to calculate QC metrics

Panpipes uses "actions" to define in which tasks to use the provided gene list.
We encode three main actions to use gene lists to describe qualities of the cells and populate the metadata (`.obs`) of the object with the newly calculated cell QC metric.
Specify the "group" name of the genes you want to use to apply a specific action to calculate the cell QC metric

If left blank, these actions will not be performed (i.e. no calculation of % of mt genes per cell will be included in the ingestion of the data)

### Supplying custom gene lists to calculate QC metrics

The human custom genelist file can be supplied by the user in two workflows to perform the three main actions:

1. **Ingest workflow**

    pipeline_ingest config file: (pipeline.yml)

    ```yaml
    custom_genes_file: resources/qc_genelist_1.0.csv
    ```

2. **Preprocess workflow**

    pipeline_preprocess config file: (pipeline.yml)

    ```yaml
    exclude_file: resources/qc_genelist_1.0.csv
    ```

*Note that we have formatted an example file containing all genes to use in both workflows, and therefore supply the same file to both workflows but users can have independent files for each of them.*

However, if the input is from mouse data then, the custom genelist file can be supplied by the user in two workflows to perform the three main actions:

1. **Ingest workflow**

    pipeline_ingest config file: (pipeline.yml)

    ```yaml
    custom_genes_file: resources/qc_gene_list_mouse.csv
    ```

2. **Preprocess workflow**

    pipeline_preprocess config file: (pipeline.yml)

    ```yaml
    exclude_file: resources/qc_gene_list_mouse.csv
### Explaining custom gene lists actions

1. **Ingest workflow** (pipeline_ingest.py)

- **calc_proportions:** calculate proportion of reads mapping to X genes over total number of reads, per cell, using [scanpy.pp.calculate_qc_metrics](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics).

    For example, for the rna modality, including a list of mitochondrial
    genes in the group `mt`, and setting

        calc_proportion: mt 

    will calculate the proportion of reads mapping to the genes whose group is "mt" and and will add `pct_counts_mt`, and `total_counts_mt`
    to the `mdata["rna"].obs` assay.

- **score_genes:** using [scanpy.tl.score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html), the action calculates the average expression of a set of genes, subtracted of the average expression of a reference set of genes. First introduced in Satija et al. Nature Biotechnology (2015).
  
  This action will generate a column **GroupName_score** and add it to the `mdata[MOD].obs` assay. For example for the rna modality, including a list of Markers for a cell type in the group 'MarkersNeutro' and setting

      score_genes: MarkersNeutro
  
  will score the cells for the list of genes provided and add a column `MarkersNeutro_score` to the `mdata["rna"].obs` assay.

2. **Preprocess workflow** (pipeline_preprocess.py)

- **exclude:** exclude these genes from the HVG selection, if they are deemed Highly Variable.

    For the exclude action, if set to `default` the workflow will look for genes whose group is set to `exclude` in the supplied qc_genelist file. Alternatively, if you are specifying your custom gene list and you want to exclude another set of genes, for example a group you call `TCR_genes`, specify this group (i.e. `exclude: TCR_genes`)

### Cell cycle actions

As described before, we also rely on a user-supplied list of genes to calculate the cell cycle phase of a cell. We believe that this choice offers the maximum flexibility to use a trusted gene-set for the calculation of this metric.
The cell cycle scoring happens in the `ingest` workflow using the `ccgenes` parameter. The cell cycle action is performed using [`scanpy.tl.score_genes_cell_cycle`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes_cell_cycle.html)

**ccgenes:**  
Setting the `ccgenes` param to `default` in the ingest workflow will calculate the phase of the cell cycle in which the cell is by using `scanpy.tl.score_genes_cell_cycle` using the file provided in panpipes/resources/cell_cicle_genes.tsv. Using this file, this action will produce at least 3 columns in the `mdata["rna"].obs` assay, namely 'S_score', 'G2M_score', 'phase'.

Users can create their own list, and need to specify the path to this new file in in the `ccgenes` param to score the cells with their custom list.

If left blank, the cellcycle score will not be calculated.

Using Custom Gene lists to plot: the Visualization workflow
---------------

Users may also supply custom gene lists to plot markers using standard visualizations, such as UMAPs or dotplots of gene expressions.
We have designated an entire workflow just for this purpose. The Visualization workflow accepts custom gene files in the [same 3-column format as defined above](#custom-gene-lists).

These files can be specified in the `viz` configuration file as follows:

pipeline_vis config file: (pipeline.yml)

```yaml
# the full list will be plotted in dot plots and matrix plots, one plot per group
full:
 - long_file1.csv
 - long_file2.csv
# the shorter list will be plotted on umaps as well as other plot types, one plot per group
minimal:
 - short_file1.csv

```

Generally in the visualization pipeline all gene groups in the input are plotted. In heatmaps and dotplots, one dotplot per group is plotted. For UMAPs, one plot per gene is
plotted, and a new file is saved per group.

## Final notes

Be deliberate and informative with the choice of group names for any gene set use, since the `.obs` column generated as output will be named based on the group of the gene list input file.
The columns added to the `MuData` object will be used for filtering in the `preprocess` workflow (see [filtering instructions](./filter_dict_instructions.md) ).

If the mitochondrial genes are in group "mt" as in the example given in the resource file, then the column generated with the **calc_proportions** action and containing the percentage of MT genes will be named **"pct_counts_mt"**.

So specifying *pct_counts_mito* instead of *pct_counts_mt* will not filter the mudata based on mitochondrial %, because the workflow can't find the supplied column in the data.
