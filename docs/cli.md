# lanctools CLI

`lanctools` provides a command-line interface with two utilities for file
conversions.

---

## Basic Usage

To get help, use the `--help` flag:

```bash
lanctools --help
lanctools merge --help
```

---

## merge

The `merge` command is used to combine a list of `.lanc` files into a single
`.lanc` file.

### Example

```bash
lanctools merge --input chr1.lanc --input chr2.lanc --input chr3.lanc --output chr1_3.lanc
```

### Arguments

| Option      | Argument | Type | Description |
| --- | --- | --- | --- |
| `--input` | TEXT | optional | A local ancestry .lanc file to be merged. This option should be repeated to specify multiple files. |
| `--input-list` | TEXT | optional | File containing local ancestry .lanc files to be merged, one per line |
| `--output` | TEXT | required | The output .lanc file |


!!! warning
    
    Either `--input` or `--input-list` may be provided, not both.

---

## convert

The `convert` command is used to convert output from a local ancestry
inference tool into the `.lanc` format. Currently supported formats
are FLARE and RFMix.

### Example

```bash
lanctools convert --input chr1.anc.vcf.gz --plink chr1 --format FLARE --output chr1.lanc
```

### Arguments


| Option      | Argument | Type | Description |
| --- | --- | --- | --- |
| `--input` | TEXT | required | The input local ancestry file |
| `--plink` | TEXT | required | The corresponding plink2 file prefix |
| `--format` | TEXT | required | The input file format ("RFMix" or "FLARE") |
| `--output` | TEXT | required | The output .lanc file |


!!! info

    Since not all positions in the plink2 files may exist in the FLARE input files, it is necessary to interpolate between positions. We use linear interpolation. This means that if e.g. the `.vcf.gz` file reports ancestry A for a given haplotype at chr1:100 and ancestry B at chr1:200, we assign ancestry A to all variants in chr1:100-150 and B to all variants from in chr1:150-200.
