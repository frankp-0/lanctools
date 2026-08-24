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

from contextlib import ExitStack
from pathlib import Path

import numba as nb
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from pandas import DataFrame
from pgenlib import PgenReader, PvarReader

from ._cpp import read_flare, read_rfmix

### ─────────────────────────────────────────────────────────────
### Functions
### ─────────────────────────────────────────────────────────────


def _parse_lanc_line(
    line: str,
) -> tuple[NDArray[np.uint8], NDArray[np.uint8], NDArray[np.uint32]]:
    """Parse a single line of .lanc file into a tuple of ancestries and breakpoints"""
    fields = line.strip().split()
    if not fields:
        raise ValueError("A .lanc sample line must contain at least one tract")

    breakpoints, left_haps, right_haps = [], [], []
    previous_breakpoint = -1
    for field in fields:
        parts = field.split(":")
        if len(parts) != 2 or len(parts[1]) != 2:
            raise ValueError(f"Invalid .lanc tract: {field!r}")

        try:
            breakpoint = int(parts[0])
            left_hap = int(parts[1][0])
            right_hap = int(parts[1][1])
        except ValueError as exc:
            raise ValueError(f"Invalid .lanc tract: {field!r}") from exc

        if breakpoint <= previous_breakpoint:
            raise ValueError("Breakpoints in a .lanc sample line must be increasing")
        if not 0 <= breakpoint <= np.iinfo(np.uint32).max:
            raise ValueError(f"Breakpoint is outside uint32 range: {breakpoint}")
        if not 0 <= left_hap <= np.iinfo(np.uint8).max:
            raise ValueError(f"Ancestry value is outside uint8 range: {left_hap}")
        if not 0 <= right_hap <= np.iinfo(np.uint8).max:
            raise ValueError(f"Ancestry value is outside uint8 range: {right_hap}")

        breakpoints.append(breakpoint)
        left_haps.append(left_hap)
        right_haps.append(right_hap)
        previous_breakpoint = breakpoint
    return (
        np.array(left_haps, np.uint8),
        np.array(right_haps, np.uint8),
        np.array(breakpoints, np.uint32),
    )


def _parse_lanc_header(line: str) -> tuple[int, int]:
    fields = line.strip().split()
    if len(fields) != 2:
        raise ValueError("A .lanc header must contain variant and sample counts")
    try:
        n_variants, n_samples = (int(field) for field in fields)
    except ValueError as exc:
        raise ValueError(f"Invalid .lanc header: {line.strip()!r}") from exc
    if n_variants <= 0 or n_samples <= 0:
        raise ValueError("A .lanc header must contain positive counts")
    return n_variants, n_samples


def _read_lanc(path: str | Path) -> FlatLanc:
    """Read a .lanc file into a FlatLanc object"""
    left_haps, right_haps, breakpoints, offsets = [], [], [], [0]
    with open(path) as f:
        try:
            n_variants, n_samples = _parse_lanc_header(next(f))
        except StopIteration as exc:
            raise ValueError("A .lanc file must contain a header") from exc

        for sample_idx in range(n_samples):
            try:
                line = next(f)
            except StopIteration as exc:
                raise ValueError(
                    f".lanc header declares {n_samples} samples, found {sample_idx}"
                ) from exc
            left_hap, right_hap, end = _parse_lanc_line(line)
            left_haps.append(left_hap)
            right_haps.append(right_hap)
            breakpoints.append(end)
            offsets.append(offsets[-1] + len(end))

        if any(line.strip() for line in f):
            raise ValueError(".lanc file contains more sample lines than its header declares")

    left_haps_all = np.concatenate(left_haps)
    right_haps_all = np.concatenate(right_haps)
    breakpoints_all = np.concatenate(breakpoints)
    return FlatLanc(
        left_haps_all,
        right_haps_all,
        breakpoints_all,
        np.array(offsets, dtype=np.uint32),
        n_variants=n_variants,
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
    if not files:
        raise ValueError("At least one .lanc file is required")

    with ExitStack() as stack:
        fs = [stack.enter_context(open(file)) for file in files]
        try:
            headers = [_parse_lanc_header(next(f)) for f in fs]
        except StopIteration as exc:
            raise ValueError("A .lanc file must contain a header") from exc
        nvars = [header[0] for header in headers]
        offset_bp = np.insert(np.cumsum(nvars)[:-1], 0, 0)
        nsamps = [header[1] for header in headers]
        if len(set(nsamps)) > 1:
            raise ValueError("Files have different numbers of samples")
        nsamp = nsamps[0]
        nvar = sum(nvars)
        out_lines = []
        for sample_idx in range(nsamp):
            try:
                lines = [_parse_lanc_line(next(f)) for f in fs]
            except StopIteration as exc:
                raise ValueError(
                    f".lanc file ended before sample {sample_idx} was complete"
                ) from exc
            for file_idx, (line, declared_nvar) in enumerate(zip(lines, nvars, strict=True)):
                if int(line[2][-1]) != declared_nvar:
                    raise ValueError(
                        f"Input file {files[file_idx]!r} has a final breakpoint that "
                        "does not match its declared variant count"
                    )
            left_haps = np.concatenate([line[0] for line in lines])
            right_haps = np.concatenate([line[1] for line in lines])
            breakpoints = np.concatenate([lines[j][2] + offset_bp[j] for j in range(len(lines))])
            if np.any(np.diff(breakpoints) <= 0):
                raise ValueError("Merged breakpoints must be strictly increasing")
            linelist = [
                f"{bp}:{hap0}{hap1}"
                for bp, hap0, hap1 in zip(breakpoints, left_haps, right_haps, strict=True)
            ]
            out_lines.append(" ".join(linelist))

        if any(line.strip() for f in fs for line in f):
            raise ValueError(".lanc file contains more sample lines than its header declares")

    with open(outfile, "w") as f:
        f.write(f"{nvar} {nsamp}\n" + "\n".join(out_lines) + "\n")


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

    required_columns = {"sample", "chrom", "spos", "epos", "anc0", "anc1"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"Local ancestry input is missing columns: {missing}")
    if df.empty:
        raise ValueError("Local ancestry input contains no ancestry tracts")

    ## Read plink files
    pvar = PvarReader(bytes(plink_prefix + ".pvar", "utf8"))
    n_variants = pvar.get_variant_ct()
    if n_variants == 0:
        raise ValueError("PLINK input contains no variants")

    ## Variant plink info
    df_pvar = _get_info(pvar, np.arange(n_variants))  # variant info
    chr_order = df_pvar["CHR"].unique()
    unknown_chromosomes = set(df["chrom"]) - set(chr_order)
    if unknown_chromosomes:
        unknown = ", ".join(sorted(str(chrom) for chrom in unknown_chromosomes))
        raise ValueError(f"Local ancestry input contains unknown chromosomes: {unknown}")
    df["chrom"] = pd.Categorical(df["chrom"], categories=chr_order, ordered=True)

    ## Sample plink info
    n_skip = 0
    with open(plink_prefix + ".psam") as psam:
        for line in psam:
            if line.startswith("#IID") | line.startswith("#FID"):
                break
            n_skip += 1

    df_psam = pd.read_csv(plink_prefix + ".psam", sep="\\s+", skiprows=n_skip, dtype=str)
    if "#IID" not in df_psam:
        raise ValueError("PLINK .psam input must contain a #IID column")
    if df_psam["#IID"].duplicated().any():
        raise ValueError("PLINK .psam input contains duplicate sample IDs")
    samples = df_psam["#IID"]

    input_samples = set(df["sample"])
    plink_samples = set(samples)
    missing_samples = plink_samples - input_samples
    extra_samples = input_samples - plink_samples
    if missing_samples or extra_samples:
        details = []
        if missing_samples:
            details.append(f"missing {sorted(missing_samples)}")
        if extra_samples:
            details.append(f"unexpected {sorted(extra_samples)}")
        raise ValueError("PLINK and local ancestry samples differ: " + "; ".join(details))

    ## Filter input to ordered plink samples
    df = df[df["sample"].isin(samples)].copy()

    ## Sort df by sample, chrom, spos
    df["sample"] = pd.Categorical(df["sample"], categories=samples, ordered=True)
    df = df.sort_values(by=["sample", "chrom", "spos"]).reset_index(drop=True)  # pyright: ignore[reportCallIssue]

    ## Exclude tracts starting after or ending before pgen range
    min_pvar = int(np.min(df_pvar["BP"]))
    max_pvar = int(np.max(df_pvar["BP"]))
    tracts_mask = (df["spos"] < max_pvar) & (df["epos"] > min_pvar)
    df = df[tracts_mask]
    if df.empty:
        raise ValueError("No ancestry tracts overlap the PLINK variant range")

    ## Clip tracts positions to pgen start, end
    df.loc[df["epos"] > max_pvar, "epos"] = max_pvar
    df.loc[df["spos"] < min_pvar, "spos"] = min_pvar

    ## Get index of first pvar pos >= tract epos
    df["idx"] = -1
    for chrom, df_chr in df.groupby("chrom", sort=False, observed=True):
        pvar_chr = df_pvar[df_pvar["CHR"] == chrom]

        local_idx = np.searchsorted(
            np.asarray(pvar_chr["BP"], dtype=int),
            np.asarray(df_chr["epos"], dtype=int),
            side="right",
        )

        # Convert chromosome-local position to global df_pvar position
        offset = pvar_chr.index[0]

        df.loc[df_chr.index, "idx"] = offset + local_idx

    ## If multiple tracts have same idx, pick last one
    df = (
        df.sort_values(["sample", "idx"])  # pyright: ignore[reportCallIssue]
        .groupby(["sample", "idx"], as_index=False, observed=True)
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
        df.groupby(["sample"], observed=True)["switch"]
        .apply(lambda x: " ".join(x.astype(str)))
        .reset_index(drop=True)
    )

    ## Write output
    header = f"{len(df_pvar)} {len(df_psam)}"
    with open(output, "w") as f:
        f.write(header + "\n" + "\n".join(lines.astype(str)) + "\n")


@nb.njit(parallel=True, cache=True)
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


def _get_geno(pgen: PgenReader, indices: NDArray[np.unsignedinteger]) -> NDArray[np.int32]:
    """Query genotypes"""
    n = pgen.get_raw_sample_ct()
    v = len(indices)
    alleles = np.empty((v, 2 * n), dtype=np.int32)
    pgen.read_alleles_list(indices, alleles)
    return alleles.reshape(v, n, 2).transpose(1, 0, 2)


def _validate_indices(indices: NDArray[np.integer], n_variants: int) -> NDArray[np.integer]:
    indices = np.asarray(indices)
    if indices.ndim != 1:
        raise ValueError("Variant indices must be a one-dimensional array")
    if not np.issubdtype(indices.dtype, np.integer):
        raise TypeError("Variant indices must have an integer dtype")
    if np.issubdtype(indices.dtype, np.signedinteger) and np.any(indices < 0):
        raise IndexError("Variant indices must be non-negative")
    if np.any(indices >= n_variants):
        raise IndexError("Variant index is outside the PLINK variant range")
    return indices


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
        n_variants: int | None = None,
    ):
        arrays = (left_haps, right_haps, breakpoints, offsets)
        if any(np.asarray(array).ndim != 1 for array in arrays):
            raise ValueError("FlatLanc arrays must be one-dimensional")
        if not (len(left_haps) == len(right_haps) == len(breakpoints)):
            raise ValueError("FlatLanc tract arrays must have equal lengths")
        if len(offsets) == 0 or offsets[0] != 0 or offsets[-1] != len(breakpoints):
            raise ValueError("FlatLanc offsets must delimit the tract arrays")
        if np.any(np.diff(offsets) < 0):
            raise ValueError("FlatLanc offsets must be non-decreasing")
        if n_variants is not None and n_variants <= 0:
            raise ValueError("FlatLanc variant count must be positive")

        self.left_haps = left_haps
        self.right_haps = right_haps
        self.breakpoints = breakpoints
        self.offsets = offsets
        self.n_variants = n_variants

    def get_lanc(self, indices: NDArray[np.unsignedinteger]) -> NDArray[np.uint8]:
        """Query phased local ancestry.

        Args:
            indices: The variant indices in (0-based)

        Returns:
            An array of ancestries, shape (N, V, 2)
        """

        if self.n_variants is not None:
            indices = _validate_indices(indices, self.n_variants)
        else:
            indices = np.asarray(indices)
            if indices.ndim != 1:
                raise ValueError("Variant indices must be a one-dimensional array")
            if not np.issubdtype(indices.dtype, np.integer):
                raise TypeError("Variant indices must have an integer dtype")
            if np.issubdtype(indices.dtype, np.signedinteger) and np.any(indices < 0):
                raise IndexError("Variant indices must be non-negative")
            if np.any(indices >= self.breakpoints[-1]):
                raise IndexError("Variant index is outside the .lanc variant range")

        idx_order = np.argsort(indices)
        idx_ordered = np.ascontiguousarray(indices[idx_order])
        idx_inverse = np.argsort(idx_order)
        left, right = _get_lanc(
            self.left_haps,
            self.right_haps,
            self.breakpoints,
            self.offsets,
            idx_ordered,
        )
        return np.stack((left[:, idx_inverse], right[:, idx_inverse]), axis=-1)


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
        ancestries: list[str] | None = None,
    ):
        """Constructs a LancData from plink2 files.

        Args:
            plink_prefix (str): The prefix for a plink2 fileset.
            lanc_file (str): The path to a .lanc file.
            ancestries (Optional[list[str]): An optional list of ordered ancestry
                names corresponding to the .lanc file.
        """
        with ExitStack() as stack:
            pgen = PgenReader(bytes(plink_prefix + ".pgen", "utf8"))
            stack.callback(pgen.close)
            pvar = PvarReader(bytes(plink_prefix + ".pvar", "utf8"))
            stack.callback(pvar.close)
            lanc = _read_lanc(lanc_file)
            n_variants = pvar.get_variant_ct()
            if lanc.n_variants != n_variants:
                raise ValueError("PLINK and .lanc files have different numbers of variants")
            if lanc.offsets.shape[0] - 1 != pgen.get_raw_sample_ct():
                raise ValueError("PLINK and .lanc files have different numbers of samples")

            if ancestries is None:
                all_values = np.concatenate([lanc.left_haps, lanc.right_haps])
                ancestries = [str(i) for i in np.unique(all_values)]
            elif len(ancestries) <= int(max(lanc.left_haps.max(), lanc.right_haps.max())):
                raise ValueError("Ancestry names do not cover all values in the .lanc file")
            stack.pop_all()

        self.pgen = pgen
        self.pvar = pvar
        self.lanc = lanc
        self.ancestries = ancestries
        self.plink_prefix = plink_prefix
        self._closed = False

    def close(self) -> None:
        """Release the underlying PLINK readers.

        Calling ``close`` more than once is safe. Query methods cannot be used
        after the readers have been closed.
        """
        if not self._closed:
            self.pgen.close()
            self.pvar.close()
            self._closed = True

    def __enter__(self) -> LancData:
        if self._closed:
            raise RuntimeError("LancData is closed")
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    def _ensure_open(self) -> None:
        if self._closed:
            raise RuntimeError("LancData is closed")

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

        self._ensure_open()
        return _get_info(self.pvar, _validate_indices(indices, self.pvar.get_variant_ct()))

    def get_lanc(self, indices: NDArray[np.unsignedinteger]) -> NDArray[np.uint8]:
        """Query phased local ancestry.

        Args:
            indices: The variant indices in pvar order (0-based), shape (V,)

        Returns:
            An array of ancestries, shape (N, V, 2)
        """

        self._ensure_open()
        return self.lanc.get_lanc(_validate_indices(indices, self.pvar.get_variant_ct()))

    def get_lanc_dosage(self, indices: NDArray[np.uint32]) -> NDArray[np.int32]:
        """Query local ancestry dosage.

        Args:
            indices: An array of variant indices in pvar order (0-based), shape (V,)

        Returns:
            An array of local ancestry dosages, shape (N, V, K) (where K is the
                number of ancestries)
        """

        self._ensure_open()
        indices = _validate_indices(indices, self.pvar.get_variant_ct())
        lanc = np.asarray(self.get_lanc(indices), dtype=np.uint8)
        ancestries = np.arange(len(self.ancestries), dtype=np.uint8)
        left_haps_mask = (lanc[:, :, 0:1] == ancestries[None, None, :]).astype(np.int32)
        right_haps_mask = (lanc[:, :, 1:2] == ancestries[None, None, :]).astype(np.int32)
        return left_haps_mask + right_haps_mask

    def get_geno(self, indices: NDArray[np.uint32]) -> NDArray[np.int32]:
        """Query phased genotypes.

        Args:
            indices: An array of variant indices (0-based)

        Returns:
            An array of phased genotypes, shape (N, V, 2)
        """

        self._ensure_open()
        return _get_geno(
            self.pgen,
            _validate_indices(indices, self.pvar.get_variant_ct()),
        )

    def get_lanc_geno(self, indices: NDArray[np.unsignedinteger]) -> NDArray[np.int32]:
        """Query genotypes deconvoluted/masked by ancestry.

        Args:
            indices: An array of variant indices (0-based)

        Returns:
            An array of genotypes masked by ancestry, shape (N, V, 2)
        """
        self._ensure_open()
        indices = _validate_indices(indices, self.pvar.get_variant_ct())
        geno = np.asarray(self.get_geno(indices), dtype=np.int32)
        lanc = np.asarray(self.lanc.get_lanc(indices), dtype=np.uint8)
        ancestries = np.arange(len(self.ancestries), dtype=np.uint8)
        left_haps_mask = (lanc[:, :, 0:1] == ancestries[None, None, :]).astype(np.int32)
        right_haps_mask = (lanc[:, :, 1:2] == ancestries[None, None, :]).astype(np.int32)
        geno_masked = left_haps_mask * geno[:, :, 0:1] + right_haps_mask * geno[:, :, 1:2]
        return geno_masked
