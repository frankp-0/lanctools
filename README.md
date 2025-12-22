# lanctools

Tools for working with local ancestry data in the `.lanc` file format.
This package contains two main components:

1. A function (and corresponding CLI) for converting FLARE and RFMix files to
the .lanc format
2. A `LancData` class with methods for fast querying of local ancestry and
local ancestry-masked genotypes

## Installation

```
pip install lanctools
```

## Quickstart

To load and query local ancestry data for a set of variants::

```
import numpy as np

from lanctools import LancData

ld = LancData(
    plink_prefix="chr1",
    lanc_file="chr1.lanc",
    ancestries=["YRI", "CEU"]
)

idx_var = np.arange(100, dtype=np.uint32)

lanc = ld.get_lanc(idx_var) # (N, 100, 2): phased local ancestry
geno = ld.get_geno(idx_var) # (N, 100, 2): phased genotypes
lanc_geno = ld.get_lanc_geno(idx_var) # (N, 100, len(ancestries))
```

To convert a FLARE (or RFMix) local ancestry file to `.lanc`:

```
from lanctools import convert_to_lanc

convert_to_lanc(
    file="chr1.anc.vcf.gz",
    file_fmt="FLARE",
    plink_prefix="chr1",
    output="chr1.lanc"
)
```

To perform the above conversion using the CLI tool:

```
lanctools convert-flare --file chr1.anc.vcf.gz --plink_prefix chr1 --output chr1.lanc
```
