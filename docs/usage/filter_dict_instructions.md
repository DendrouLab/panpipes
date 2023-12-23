# Filtering mudata objects within preprocess pipeline

In `panpipes preprocess` a completely customisable filtering process is implemented, such that you can filter your data by basically any metric.

The filtering process in panpipes is sequential as it goes through the filtering dictionary in the `pipeline.yml`.

Within the yaml file the filtering dictionary has this basic structure:

```yaml
rna/prot/atac:
  obs:
    max:
    min:
    bool:
  var:
    max:
    min:
    bool:
  
```

For each modality, starting with rna, panpipes will first filter on `obs` and then `var`.
Each modality has a dictionary which must be fully customised to the columns in your the mudata.obs or var object.

**When specifying a column name, make sure it exactly matches the column name in the h5mu object.** 

You can review this by loading up your h5mu object in Python:

```python
import muon as mu; mu.read(filepath).obs.columns
```

Firstly for columns in `mdata['mod'].obs`.
Under "max" you list any columns which you want to run a maximum filter for e.g. total_counts, under "min" you list any columns you want to run a minimum filter and under "bool" list any boolean columns you want to filter on e.g. "is_doublet"
Then repeat for columns in `mdata['mod'].var`.
For example:

```yaml
rna:
  obs:
    min:  
      n_genes_by_counts: 500 
    max: 
      total_counts: 20000
      pct_counts_mt: 20
    bool:
      is_doublet: False
  var:
    min:
      n_cells_by_counts: 3
```

In this example, the cells are filtered to contain more than 500 genes, less than 20000 counts, less than 20% mitochondrial content, and where is_doublet is False. Then genes are filtered to contain at least 3 cells with >0 reads.

You can use this notation to expand filtering to any set of columns you have in any modality.
