# Python API

**lanctools** exports two classes for working with local ancestry data. [LancData](#lanctools.LancData) contains the genotype and local ancestry data for a set of plink2 `.pgen` and `.lanc` files, together with efficient methods for querying this data. [FlatLanc](#lanctools.FlatLanc) is the core data structure which stores local ancestry data in a flattened structure.

`LancData` owns open PLINK readers. Use it as a context manager, or call
`close()` when finished, especially when processing many datasets:

```python
with LancData(plink_prefix="chr1", lanc_file="chr1.lanc") as data:
    lanc = data.get_lanc(indices)
```

## ::: lanctools.LancData

    handler: python
    options:
      show_root_heading: true
      show_source: true

## ::: lanctools.FlatLanc

    handler: python
    options:
      show_root_heading: true
      show_source: true
