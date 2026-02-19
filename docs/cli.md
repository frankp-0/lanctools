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

```{bash}
lanctools merge --file chr1.lanc,chr2.lanc,chr3.lanc --outfile chr1_3.lanc
```

### Arguments

| Option      | Argument | Type | Description |
| --- | --- | --- | --- |
| `--file` | TEXT | required | A comma-separated list of input .lanc files |
| `--outfile` | TEXT | required | The output .lanc file |


---

## convert-flare

The `convert-flare` command is used to convert `.vcf.gz` files with annotations
for local ancestry (as output by FLARE) into the `.lanc` file format.

### Example

```bash
lanctools convert-flare --file chr1.anc.vcf.gz --plink-prefix chr1 --output chr1.lanc
```

### Arguments


| Option      | Argument | Type | Description |
| --- | --- | --- | --- |
| `--file` | TEXT | required | The input file, or a list of comma-separated input files |
| `--plink-prefix` | TEXT | required | The corresponding plink2 file prefix or comma-separated list of prefixes |
| `--lanc-file` | TEXT | required | The output .lanc file or comma-separated list of output files |


!!! info

    Since not all positions in the plink2 files may exist in the FLARE input
    files, it is necessary to interpolate between positions. We use linear
    interpolation. This means that if e.g. the `.vcf.gz` file reports
    ancestry A for a given haplotype at chr1:100 and ancestry B at chr1:200, we
    assign ancestry A to all variants in chr1:100-150 and B to all variants
    from in chr1:150-200.

!!! warning
    
    All arguments must have the same length. E.g. if two comma-separated FLARE
    files are provided to `--file`, two files must be provided to `--plink-prefix`
    and `--lanc-file`.

---

## convert-rfmix

The `convert-flare` command is used to convert `.msp.tsv` files output by RFMix
into the `.lanc` file format.

### Example

```bash
lanctools convert-rfmix --file chr1.msp.tsv --plink-prefix chr1 --output chr1.lanc
```

### Arguments


| Option      | Argument | Type | Description |
| --- | --- | --- | --- |
| `--file` | TEXT | required | The input file, or a list of comma-separated input files |
| `--plink-prefix` | TEXT | required | The corresponding plink2 file prefix or comma-separated list of prefixes |
| `--lanc-file` | TEXT | required | The output .lanc file or comma-separated list of output files |


!!! warning
    
    All arguments must have the same length. E.g. if two comma-separated RFMix
    files are provided to `--file`, two files must be provided to `--plink-prefix`
    and `--lanc-file`.

