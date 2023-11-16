Reference Mapping (refmap)
======

The reference mapping pipeline `panpipes refmap` implements scvi-tools reference mapping tools `scvi`, `totalvi` and `scanvi`.
You can supply an query anndata/mudata containing RNA, and the reference model `path_to_model\model.pt`. Since reference and query have to work from the same subset of genes, you can also supply the reference anndata/mudata containing RNA or RNA+PROT (for `totalvi`) so the genes selection can be unified and the resulting plots include the reference cells together with the query. 
Note that you can even use a reference scvi model created previously by `panpipes integration`.

Steps:
------

1.  In a new folder, generate config file for integration,
    `panpipes refmap config` and edit the pipeline.yml file.
2.  Run complete refmap pipeline with `panpipes refmap make full`


## Expected structure of MuData object
The ideal way to run `panpipes refmap` is to use the output mudata file from `panpipes preprocess`, as this will make sure the MuData object has correctly names layers and slots. 

## Formatting query data 

The reference mapping models used in `panpipes` require that the query is formatted in a special way to match some of the features in the reference. It is for example common practice to share the list of the highly variable genes (or the equivalent set of features) that was used to build the reference, alongside the model itself. We had already formatted the data for this tutorial, but if you're interested in using your own query you should pay special attention to this.

Since each reference model has likely a different structure, at the moment we don't provide standardised code to format the reference and query. 

For example, the reference and query data we used for the paper looks like this:

```
> reference

AnnData object with n_obs × n_vars = 152094 × 4000
    obs: 'nCount_ADT', 'nFeature_ADT', 'nCount_RNA', 'nFeature_RNA', 'orig.ident', 'lane', 'donor', 'time', 'celltype.l1', 'celltype.l2', 'celltype.l3', 'Phase', 'nCount_SCT', 'nFeature_SCT', 'X_index', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'Protein log library size', 'Number proteins detected', 'RNA log library size', '_scvi_labels', '_scvi_batch'
    var: 'mt', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches'
    uns: 'log1p', 'hvg', '_scvi_uuid', '_scvi_manager_uuid'
    obsm: 'protein_counts', 'X_totalvi'
    layers: 'counts'

> query

AnnData object with n_obs × n_vars = 57669 × 4000
    obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'RNA_snn_res.0.4', 'seurat_clusters', 'set', 'Resp', 'disease', 'subj_code', 'covidpt_orhealth', 'mito', 'ncount', 'nfeat', 'bust_21', 'og_clust', 'severmod_other', 'og_clusts', 'nCount_ADT', 'nFeature_ADT', 'UMAP1', 'UMAP2', 'final_clust', 'final_clust_v2', 'new_pt_id', 'Resp_og', 'final_clust_withnum', 'final_clust_review', 'Age', 'Gender', 'Gender_num', 'celltype.l2', '_scvi_labels', '_scvi_batch'
    uns: 'log1p', '_scvi_uuid', '_scvi_manager_uuid'
    obsm: 'pro_exp', 'protein_counts'
    layers: 'counts'

```

In this case, we had to make sure to format the query so that:
- it contains the same 4000 HVG as the reference
- it contains the protein counts in the `.obsm['protein_counts']` layer
- the protein names matched the names of the protein reference (and when not matching we padded the array with 0s to allow imputation of missing features)
- the query has a column matching the reference's celltype labels we want to transfer,  filled with `Unknown` 

```
query.obs["celltype.l2"]

AAACCCACACCAGCGT-1    Unknown
AAACCCACATCTCAAG-1    Unknown
AAACGAAAGACCTGGA-1    Unknown
AAACGCTCAGTGGGTA-1    Unknown
AAACGCTGTAGCTTGT-1    Unknown
                       ...   
TTTGGTTTCCAATCTT-1    Unknown
TTTGTTGAGACGTCCC-1    Unknown
TTTGTTGAGGTATTGA-1    Unknown
TTTGTTGCAGGCGATA-1    Unknown
TTTGTTGTCTTCTGTA-1    Unknown
Name: celltype.l2, Length: 57669, dtype: category
Categories (1, object): ['Unknown']
```

