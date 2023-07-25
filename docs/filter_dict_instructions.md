Filtering
============

The filtering process in panpipes is sequential as it goes through the filtering dictionary.

For each modality, starting with rna, it will first filter on obs and then vars.
Each modality has a dictionary which must be fully customised to the
columns in your the mudata.obs or var object. This is turn will be defined the names of the "groups" in your input gene list files. For more information review: :doc:`gene_list_format`


**When specifying a column name, make sure it exactly matches the column name in the h5mu object.** 
You can review this by loaing up your h5mu object in python: 
`python; import muon as mu; mu.read(filepath).obs.columns`

Basic dict:

```
mod:
  obs:
    max:
    min:
    bool:
  var:
    max:
    min:
    bool:
  
```

Firstly for columns in mdata['mod'].obs.
Under "max" you list any columns which you want to run a maximum filter for e.g. total_counts, under "min" you list any columns you want to run a minimum filter and under "bool" list any boolean columns you want to filter on e.g. "is_doublet"

Then repeat for columns in mdata['mod'].var

Full example for rna

```
rna:
  obs:
    min:  
      n_genes_by_counts: 500 
    max: 
      total_counts: 20000
      pct_counts_mt: 20
    bool:
      is doublet: False
  var:
    min:
      n_cells_by_counts: 3
```

In this example, the cells are filtered to contain more than 500 genes, less than 20000 counts, less than 20% mitochondrial content, and where is_doublet is False. Then genes are filtered to contain at least 3 cells with >0 reads.

You can use this notation to expand filtering to any set of columns you have in any modality. 


