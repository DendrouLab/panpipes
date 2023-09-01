Preprocessing
=============


## Pipeline steps

The preprocess pipeline filters the data as defined in the [filtering dictionary](../usage/filter_dict_instructions.md) section of the `pipeline.yml`. The data can also been downsampled.
Then each modality is normalised and scaled. For the RNA this is normalising counts per cell with [scanpy.pp.normalize_total](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html) and optionally, regressing and scaling the data using scanpy functions. Highly variable genes (HVGs) are also calculated, and a PCA performed on those highly variable genes. There is an option to exclude specific genes from the HVGs e.g. HLA genes or BCR/TCR genes. These are specified in the same way as all [gene lists](../usage/gene_list_format). In the example below, the "group" in the gene list file is "exclude".
```
hvg:
  exclude_file: resources/qc_genelist_1.0.csv
  exclude: "exclude"
```

For Protein assay, the data are normalised either by centralised-log-ratio or by dsb as described in the muon documentation [here](https://muon.readthedocs.io/en/latest/omics/citeseq.html). There is additional panpipes functionality to trim dsb outliers as discussed on the dsb [github page](https://github.com/niaid/dsb/issues/9)


For the ATAC assay ....


## Steps to run:


1. In a new folder, generate config file for integration,
   ``panpipes preprocess config``
2. edit the pipeline.yml file

   -  The filtering options are dynamic depending on your qc_mm inputs. This is described [here](../usage/filter_dict_instructions.md) 
   -  There are lots of options for normalisation explained in the
      pipeline.yml

3. Run complete preprocess pipeline with
   ``panpipes preprocess make full``

The h5mu outputted from ``preprocess`` is filtered and normalised, and
for rna highly variable genes are computed.


## Expected structure of MuData object
The ideal way to run `panpipes preprocess` is to use the output mudata file from `panpipes qc_mm`, as this will make sure the MuData object has correctly names layers and slots. 

The bare minimum MuData object required is raw data in the X slot of each modality and a sample_id column the .obs slot of each of each modality, and the common (outer) obs.
