Gene list formats
=================


In the ingestion and vis pipelines, the user can provide custom genelists
in order to compute gene list scores and to visualise


All <sup>[1](#footnote1)</sup> gene lists provided to the pipeline should be in a 3 columns format, where the column headers are "mod" (modality: "rna", "prot", or "atac"), feature and group. The group column is used to distinguish different gene groups.

mod  | feature | group
-----| --------| -------
rna  | MT-ND1  | mt
rna  | MT-ND2  | mt

The gene lists are not pre-determined within panpipes in order to maximise
flexibility as all organisms will require separate lists, but there are example lists provided on our [github page](https://github.com/DendrouLab/panpipes/tree/main/panpipes/resources)



Ingestion gene list inputs
-----------------------
In the ingestion pipeline, panpipes can score gene lists based on the percent gene content (as is typically done for mitochondrial content), or score genes based on the average expression of a set of genes subtracted with the average expression of a reference set of genes.

This is defined in the following section of the ingestion pipeline.yml:

    custom_genes_file: path/to/resources/qc_genelist_1.0.csv
    calc_proportions: hb,mt,rp
    score_genes: MarkersNeutro

Example custom_genes_file:
[resources/qc_genelist_1.0.csv](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/qc_genelist_1.0.csv).

Calc_proportions: compute
what percentage of the gene counts in each barcode are associated to the
group of genes from the provided gene list. The genes are scored for each modality using 
[scanpy.pp.calculate_qc_metrics](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics).
For example, for the rna modality, including a list of mitochondiral
genes in the group `mt`, will add `pct_counts_mt`, and `total_counts_mt`
to the `mdata["rna"].obs` assay.

Score genes: Uses the same approach as calc_proportions but instead scores a set of genes using
[scanpy.tl.score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html).

Visualisation pipeline input genelists:
--------------------

For the visualisation pipeline

pipeline.yml excerpt:

    # the full list will be plotted in dot plots and matrix plots, one plot per group
    full:
     - long_file1.csv
     - long_file2.csv
    # the shorter list will be plotted on umaps as well as other plot types, one plot per group
    minimal:
     - short_file1.csv

These require the same format as above.
Generally in the visualisation pipeline all gene groups in the input are plotted. In heatmaps and dot
plots, one dotplot per group is plotted. For umaps, one plot per gene is
plotted, and a new file is saved per group.


##### Footnotes
<a name="footnote1">1</a>:  the one exception to this is how cellcycle genes are included,
this will change in a future version of panpipes
[resources/cell_cycle_genes.csv](https://github.com/DendrouLab/panpipes/blob/main/panpipes/resources/cell_cycle_genes.tsv)