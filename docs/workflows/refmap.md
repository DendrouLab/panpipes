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