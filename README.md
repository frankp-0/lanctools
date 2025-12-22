# lanctools

Tools for working with phased local ancestry data stored in the `.lanc` file format, as defined by Admix-kit [Hou et al., 2024].

`lanctools` is designed to provide **fast local ancestry queries** and convenient conversion from external formats (e.g., FLARE [Browning et al., 2023] and RFMix [Maples et al., 2013]). It focuses on efficient access to `.lanc` data and is **not** intended to replace the full functionality of Admix-kit.

## Features

- Efficient random access to phased local ancestry data
- Local ancestry-masked genotype queries
- Conversion from FLARE and RFMix output to `.lanc` format
- Python API and command-line interface (CLI)

## Installation

```bash
pip install lanctools
```

## Quickstart

### Querying Local Ancestry Data

```python
import numpy as np
from lanctools import LancData

ld = LancData(
    plink_prefix="chr1",
    lanc_file="chr1.lanc",
    ancestries=["YRI", "CEU"]
)

idx_var = np.arange(100, dtype=np.uint32)

# Get phased local ancestry: shape (N, 100, 2)
lanc = ld.get_lanc(idx_var)

# Get phased genotypes: shape (N, 100, 2)
geno = ld.get_geno(idx_var)

# Get ancestry-masked genotypes: shape (N, 100, len(ancestries))
lanc_geno = ld.get_lanc_geno(idx_var)

```

### Converting FLARE or RFMix Files to .lanc

```python
from lanctools import convert_to_lanc

convert_to_lanc(
    file="chr1.anc.vcf.gz",
    file_fmt="FLARE",
    plink_prefix="chr1",
    output="chr1.lanc"
)
```

### Command-Line Interface

```bash
lanctools convert-flare --file chr1.anc.vcf.gz --plink_prefix chr1 --output chr1.lanc
```

## References

- Hou, K. et al. Admix-kit: an integrated toolkit and pipeline for genetic analyses of admixed populations. Bioinformatics 40, btae148 (2024). [paper](https://doi-org.libproxy.lib.unc.edu/10.1093/bioinformatics/btae148) [software](https://github.com/KangchengHou/admix-kit)
- Browning, S. R., Waples, R. K. & Browning, B. L. Fast, accurate local ancestry inference with FLARE. Am J Hum Genet 110, 326–335 (2023). [paper](https://doi.org/10.1016/j.ajhg.2022.12.010) [software](https://github.com/browning-lab/flare)
- Maples, B. K., Gravel, S., Kenny, E. E. & Bustamante, C. D. RFMix: A Discriminative Modeling Approach for Rapid and Robust Local-Ancestry Inference. Am J Hum Genet 93, 278–288 (2013). [paper](<https://doi.org/10.1016/j.ajhg.2013.06.020>) [software](https://github.com/slowkoni/rfmix)
