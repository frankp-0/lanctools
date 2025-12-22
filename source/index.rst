.. lanctools documentation master file, created by
   sphinx-quickstart on Sat Dec 20 23:18:55 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

lanctools
=========

Tools for working with local ancestry data in the `.lanc` file format.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   lanctools

Quickstart
----------

To load and query local ancestry data for a set of variants::

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
