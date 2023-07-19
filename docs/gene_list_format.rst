Gene list formats
=================

In the qc_mm and vis pipelines, the user can provide custom genelists in
order to compute gene list scores and to visualise

All [^1] gene lists provided to the pipeline should be in a 3 columns
format, where the column headers are mod, feature and group. The gorup
column is used to distinguish different gene groups.

=== ======= =====
mod feature group
=== ======= =====
rna MT-ND1  mt
rna MT-ND2  mt
=== ======= =====

The genes are not pre-determined within panpipes in order to maximise
flexibility as all organisms will require separate lists.

[^1] the one exception to this is how cellcycle genes are included, this
will change in a future version of panpipes
`resources/cell_cycle_genes.csv <https://github.com/DendrouLab/panpipes/blob/master/resources/cell_cycle_genes.csv>`__

QC_mm gene list inputs
----------------------

qc_mm pipeline.yml excerpt:

::

   custom_genes_file: path/to/resources/qc_genelist_1.0.csv
   calc_proportions: hb,mt,rp
   score_genes: MarkersNeutro

Provided custom_genes_file:
`resources/qc_genelist_1.0.csv <https://github.com/DendrouLab/panpipes/blob/master/resources/qc_genelist_1.0.csv>`__

Calc proportions: will take which ever groups are specified and compute
what percentage of the gene counts in each barcode are associated to the
list. This is done by listibng them as qc_vars in
`scanpy.pp.calculate_qc_metrics <https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics>`__.
For example, for the rna modality, including a list of mitochondiral
genes in the group “mt”, will add pct_counts_mt, and total_counts_mt to
the mdata["rna"].obs assay.

Score genes: score a set of genes using
`scanpy.tl.score_genes <https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html>`__

Vis input genelists:
--------------------

For the visualisation pipeline

pipeline.yml excerpt:

::

   # the full list will be plotted in dot plots and matrix plots, one plot per group
   full:
    - long_file1.csv
    - long_file2.csv
   # the shorter list will be plotted on umaps as well as other plot types, one plot per group
   minimal:
    - short_file1.csv

These require the same format as above, except that you do not specify
which groups are plot, all groups are plotted. In heatmaps and dot
plots, one dotplot per group is plotted. For umaps, one plot per gene is
plotted, and a new file is saved per group.
