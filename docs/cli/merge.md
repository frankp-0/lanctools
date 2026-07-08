# Merging .lanc Files

The `merge` command is used to combine a list of `.lanc` files into a single
`.lanc` file.

## Example

```bash
lanctools merge --input chr1.lanc --input chr2.lanc --input chr3.lanc --output chr1_3.lanc
```

## Arguments

| Option      | Argument | Type | Description |
| --- | --- | --- | --- |
| `--input` | TEXT | optional | A local ancestry .lanc file to be merged. This option should be repeated to specify multiple files. |
| `--input-list` | TEXT | optional | File listing local ancestry .lanc files to be merged, one per line |
| `--output` | TEXT | required | The output .lanc file |

!!! warning

    Either `--input` or `--input-list` may be provided, not both.
