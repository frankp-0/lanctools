# File Conversions

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

!!! note

    For `--format RFMix`, the `--input` file should be RFMix's `.msp.tsv` output file, which contains the most likely assignment of subpopulations per CRF point. For `--format FLARE`, the `--input` file should be the `.anc.vcf.gz` file output by FLARE.

!!! info

    Since not all positions in the plink2 files may exist in the FLARE input files, it is necessary to interpolate between positions. We use linear interpolation. This means that if e.g. the `.vcf.gz` file reports ancestry A for a given haplotype at chr1:100 and ancestry B at chr1:200, we assign ancestry A to all variants in chr1:100-150 and B to all variants from in chr1:150-200.
