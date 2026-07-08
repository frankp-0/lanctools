# MIT License
# Copyright (c) 2025 Franklin Ockerman
# See LICENSE file for full license text

"""Core module for lanctools.

This module provides:

- Conversion from other local ancestry formats to .lanc with `convert_to_lanc`
- Merging multiple .lanc  files with `merge_lanc`
- The `FlatLanc` class, which represents local ancestry data in a flattened
  structure for fast querying using `_get_lanc`
- The `LancData` class, which provides an interface for querying local
  ancestry, genotypes, and ancestry-deconvoluted genotypes
"""

from __future__ import annotations
from pathlib import Path
from pgenlib import PgenReader, PvarReader
import numpy as np
from numpy.typing import NDArray
import numba as nb
import pandas as pd
from pandas import DataFrame
from typing import Optional
from ._cpp import read_rfmix, read_flare


### ─────────────────────────────────────────────────────────────
### Functions
### ─────────────────────────────────────────────────────────────


def _parse_lanc_line(
    line: str,
) -> tuple[NDArray[np.uint8], NDArray[np.uint8], NDArray[np.uint32]]:
    """Parse a single line of .lanc file into a tuple of ancestries and breakpoints"""
    fields = line.strip().split()
    breakpoints, left_haps, right_haps = [], [], []
    for field in fields:
        breakpoint, hap_pair = field.split(":")
        breakpoints.append(int(breakpoint))
        left_haps.append(int(hap_pair[0]))
        right_haps.append(int(hap_pair[1]))
    return (
        np.array(left_haps, np.uint8),
        np.array(right_haps, np.uint8),
        np.array(breakpoints, np.uint32),
    )


def _read_lanc(path: str | Path) -> FlatLanc:
    """Read a .lanc file into a FlatLanc object"""
    left_haps, right_haps, breakpoints, offsets = [], [], [], [0]
    with open(path, "r") as f:
        next(f)
        for line in f:
            left_hap, right_hap, end = _parse_lanc_line(line)
            left_haps.append(left_hap)
            right_haps.append(right_hap)
            breakpoints.append(end)
            offsets.append(offsets[-1] + len(end))

    left_haps_all = np.concatenate(left_haps)
    right_haps_all = np.concatenate(right_haps)
    breakpoints_all = np.concatenate(breakpoints)
    return FlatLanc(
        left_haps_all,
        right_haps_all,
        breakpoints_all,
        np.array(offsets, dtype=np.uint32),
    )


def _get_info(pvar: PvarReader, indices: NDArray[np.unsignedinteger]) -> DataFrame:
    chrom = [pvar.get_variant_chrom(i).decode("utf8") for i in indices]
    pos = [pvar.get_variant_pos(i) for i in indices]
    ref = [pvar.get_allele_code(i, 0).decode("utf8") for i in indices]
    alt = [pvar.get_allele_code(i, 1).decode("utf8") for i in indices]
    id = [pvar.get_variant_id(i).decode("utf8") for i in indices]
    df = DataFrame({"CHR": chrom, "BP": pos, "REF": ref, "ALT": alt, "ID": id})
    df["BP"] = df["BP"].astype("uint32")
    return df


def merge_lanc(files: list[str], outfile: str):
    fs = [open(f, "r") for f in files]
    line1s = [next(f).strip().split() for f in fs]
    nvars = [int(l[0]) for l in line1s]
    offset_bp = np.cumsum(nvars)[:-1]
    offset_bp = np.insert(offset_bp, 0, 0)
    nsamps = [int(l[1]) for l in line1s]
    if len(set(nsamps)) > 1:
        raise ValueError("Files have different numbers of samples")
    nsamp = int(nsamps[0])
    nvar = sum(nvars)
    out_lines = []
    for i in range(nsamp):
        lines = [_parse_lanc_line(next(f)) for f in fs]
        left_haps = np.concatenate([l[0] for l in lines])
        right_haps = np.concatenate([l[1] for l in lines])
        breakpoints = np.concatenate(
            [lines[j][2] + offset_bp[j] for j in range(len(lines))]
        )
        linelist = [
            f"{bp}:{hap0}{hap1}"
            for bp, hap0, hap1 in zip(breakpoints, left_haps, right_haps)
        ]
        line = " ".join(linelist)
        out_lines.append(line)

    ## Write output
    header = f"{nvar} {nsamp}"
    with open(outfile, "w") as f:
        f.write(header + "\n" + "\n".join(out_lines) + "\n")


def convert_to_lanc(file: str, file_fmt: str, plink_prefix: str, output: str):
    """Convert local ancestry files to .lanc format

    This function currently only supports FLARE and RFMix input.

    Args:
        file: The local ancestry file
        file_fmt: Input local ancestry format, either "FLARE" or "RFMix"
        plink_prefix: The prefix for a plink2 fileset corresonding to file
        output: The output file where the result is written
    """

    ## Read input local ancestry file to pandas DataFrame
    if file_fmt == "FLARE":
        df = pd.DataFrame(read_flare(file))
    elif file_fmt == "RFMix":
        df = pd.DataFrame(read_rfmix(file))
    else:
        raise ValueError("Please specify either `FLARE` or `RFMix` input")

    ## Read plink files
    pvar = PvarReader(bytes(plink_prefix + ".pvar", "utf8"))
    n_variants = pvar.get_variant_ct()

    ## Variant plink info
    df_pvar = _get_info(pvar, np.arange(n_variants))  # variant info

    ## Sample plink info
    n_skip = 0
    with open(plink_prefix + ".psam") as psam:
        for line in psam:
            if line.startswith("#IID") | line.startswith("#FID"):
                break
            n_skip += 1

    df_psam = pd.read_csv(
        plink_prefix + ".psam", sep="\\s+", skiprows=n_skip, dtype=str
    )
    samples = df_psam["#IID"]

    if not samples.isin(df["sample"]).all():
        raise ValueError("Not all pgen samples exist in local ancestry input")

    ## Filter input to ordered plink samples
    df = df[df["sample"].isin(samples)].copy()

    ## Sort df by sample, chrom, spos
    df["sample"] = pd.Categorical(df["sample"], categories=samples, ordered=True)
    df = df.sort_values(by=["sample", "chrom", "spos"]).reset_index(drop=True)

    ## Exclude tracts starting after or ending before pgen range
    min_pvar = np.min(df_pvar["BP"])
    max_pvar = np.max(df_pvar["BP"])
    tracts_mask = (df["spos"] < max_pvar) & (df["epos"] > min_pvar)
    df = df[tracts_mask]

    ## Clip tracts positions to pgen start, end
    df["epos"] = df["epos"].clip(upper=max_pvar)
    df["spos"] = df["spos"].clip(lower=min_pvar)

    ## Get index of first pvar pos >= tract epos
    df["idx"] = np.searchsorted(df_pvar["BP"].values, df["epos"].values, side="right")

    ## If multiple tracts have same idx, pick last one
    df = (
        df.sort_values(["sample", "chrom", "idx"])
        .groupby(["sample", "chrom", "idx"], as_index=False, observed=True)
        .tail(1)  # last row per group
    )

    ## Set ending idx of last tract to extend to end of pvar
    idxmax_rows = df.groupby("sample", observed=True)["idx"].idxmax()
    df.loc[idxmax_rows, "idx"] = len(df_pvar)

    ## Get .lanc file lines
    df["switch"] = (
        df["idx"]
        .astype(str)
        .str.cat(df["anc0"].astype(str), sep=":")
        .str.cat(df["anc1"].astype(str))
    )
    lines = (
        df.groupby(["sample", "chrom"], observed=True)["switch"]
        .apply(lambda x: " ".join(x.astype(str)))
        .reset_index(drop=True)
    )

    ## Write output
    header = f"{len(df_pvar)} {len(df_psam)}"
    with open(output, "w") as f:
        f.write(header + "\n" + "\n".join(lines.astype(str)) + "\n")


@nb.njit(parallel=True)
def _get_lanc(
    left_haps: NDArray[np.uint8],
    right_haps: NDArray[np.uint8],
    breakpoints: NDArray[np.uint32],
    offsets: NDArray[np.uint32],
    indices: NDArray[np.unsignedinteger],
) -> tuple[NDArray[np.uint8], NDArray[np.uint8]]:
    """Query local ancestry"""
    n_samples = offsets.shape[0] - 1
    n_variants = indices.shape[0]
    left_out = np.empty((n_samples, n_variants), dtype=np.uint8)
    right_out = np.empty((n_samples, n_variants), dtype=np.uint8)

    for i in range(1, n_variants):
        if indices[i] < indices[i - 1]:
            raise ValueError("indices must be sorted ascending")

    for i in nb.prange(n_samples):
        start = offsets[i]
        end = offsets[i + 1]
        end_i = breakpoints[start:end]
        left_i = left_haps[start:end]
        right_i = right_haps[start:end]

        j = 0
        end_len = end_i.shape[0]
        for q in range(n_variants):
            idx = indices[q]
            while j < end_len and idx >= end_i[j]:
                j += 1
            left_out[i, q] = left_i[j]
            right_out[i, q] = right_i[j]
    return left_out, right_out


def _get_geno(
    pgen: PgenReader, indices: NDArray[np.unsignedinteger]
) -> NDArray[np.int32]:
    """Query genotypes"""
    n = pgen.get_raw_sample_ct()
    v = len(indices)
    alleles = np.empty((v, 2 * n), dtype=np.int32)
    pgen.read_alleles_list(indices, alleles)
    return alleles.reshape(v, n, 2).transpose(1, 0, 2)


### ─────────────────────────────────────────────────────────────
### Data structures
### ─────────────────────────────────────────────────────────────


class FlatLanc:
    """Stores .lanc file ancestry data in a flattened structure for fast querying.

    Attributes:
        right_haps (NDArray[uint8]): Concatenated right haplotypes for all samples, shape (H,)
        left_haps (NDArray[uint8]): Concatenated left haplotypes for all samples, shape (H,)
        breakpoints (NDArray[uint32]): Concatenated breakpoints for all samples, shape (H,)
        offsets (NDArray[uint32]): Cumulative end indices separating samples, shape (N,)
    """

    def __init__(
        self,
        left_haps: NDArray[np.uint8],
        right_haps: NDArray[np.uint8],
        breakpoints: NDArray[np.uint32],
        offsets: NDArray[np.uint32],
    ):
        self.left_haps = left_haps
        self.right_haps = right_haps
        self.breakpoints = breakpoints
        self.offsets = offsets


class LancData:
    """The genotype and local ancestry data for a single chromosome/dataset.

    Attributes:
        pgen (PgenReader): A pgenlib PgenReader object.
        pvar (PvarReader): A pgenlib PVarReader object.
        lanc (FlatLanc): A FlatLanc object with local ancestry data.
        ancestries (list[str]): An ordered list of ancestry names. The integer codes in
            the .lanc file and `self.lanc` correspond to indices in this list (e.g.
            0 -> ancestries[0]).
        plink_prefix (str): The prefix for the corresponding plink2 fileset.
    """

    def __init__(
        self,
        plink_prefix: str,
        lanc_file: str,
        ancestries: Optional[list[str]] = None,
    ):
        """Constructs a LancData from plink2 files.

        Args:
            plink_prefix (str): The prefix for a plink2 fileset.
            lanc_file (str): The path to a .lanc file.
            ancestries (Optional[list[str]): An optional list of ordered ancestry names corresponding to the .lanc file.
        """
        pgen = PgenReader(bytes(plink_prefix + ".pgen", "utf8"))
        pvar = PvarReader(bytes(plink_prefix + ".pvar", "utf8"))
        lanc = _read_lanc(lanc_file)

        if ancestries is None:
            all_values = np.concatenate([lanc.left_haps, lanc.right_haps])
            ancestries = [str(i) for i in np.unique(all_values)]

        self.pgen = pgen
        self.pvar = pvar
        self.lanc = lanc
        self.ancestries = ancestries
        self.plink_prefix = plink_prefix

    def get_info(self, indices: NDArray[np.uint32]) -> DataFrame:
        """Query info for a set of variants.

        Args:
            indices: The variant indices in pvar order (0-based), shape (V,)

        Returns:
            pandas.DataFrame: One row per variant with the following columns:

                - CHR (str): Chromosome name. \n
                - BP (int): 1-based genomic position. \n
                - REF (str): Reference allele. \n
                - ALT (str): Alternate allele. \n
                - ID (str): Variant identifier. \n
        """

        return _get_info(self.pvar, indices)

    def get_lanc(self, indices: NDArray[np.unsignedinteger]) -> NDArray[np.uint8]:
        """Query phased local ancestry.

        Args:
            indices: The variant indices in pvar order (0-based), shape (V,)

        Returns:
            An array of ancestries, shape (N, V, 2)
        """

        left, right = _get_lanc(
            self.lanc.left_haps,
            self.lanc.right_haps,
            self.lanc.breakpoints,
            self.lanc.offsets,
            indices,
        )
        return np.stack((left, right), axis=-1)

    def get_lanc_dosage(self, indices: NDArray[np.uint32]) -> NDArray[np.uint8]:
        """Query local ancestry dosage.

        Args:
            indices: An array of variant indices in pvar order (0-based), shape (V,)

        Returns:
            An array of local ancestry dosages, shape (N, V, K) (where K is the
                number of ancestries)
        """

        lanc = np.asarray(self.get_lanc(indices), dtype=np.uint8)
        ancestries = np.arange(len(self.ancestries), dtype=np.uint8)
        left_haps_mask = (lanc[:, :, 0:1] == ancestries[None, None, :]).astype(np.int32)
        right_haps_mask = (lanc[:, :, 1:2] == ancestries[None, None, :]).astype(
            np.int32
        )
        return left_haps_mask + right_haps_mask

    def get_geno(self, indices: NDArray[np.uint32]) -> NDArray[np.int32]:
        """Query phased genotypes.

        Args:
            indices: An array of variant indices in pvar order (0-based), shape (V,)

        Returns:
            An array of phased genotypes, shape (N, V, 2)
        """

        return _get_geno(self.pgen, indices)

    def get_lanc_geno(self, indices: NDArray[np.unsignedinteger]) -> NDArray[np.int32]:
        """Query genotypes deconvoluted/masked by ancestry.

        Args:
            indices: An array of variant indices in pvar order (0-based), shape (V,)

        Returns:
            An array of genotypes masked by ancestry, shape (N, V, 2)
        """
        geno = np.asarray(self.get_geno(indices), dtype=np.int32)
        lanc = np.asarray(self.get_lanc(indices), dtype=np.uint8)
        ancestries = np.arange(len(self.ancestries), dtype=np.uint8)
        left_haps_mask = (lanc[:, :, 0:1] == ancestries[None, None, :]).astype(np.int32)
        right_haps_mask = (lanc[:, :, 1:2] == ancestries[None, None, :]).astype(
            np.int32
        )
        geno_masked = (
            left_haps_mask * geno[:, :, 0:1] + right_haps_mask * geno[:, :, 1:2]
        )
        return geno_masked
