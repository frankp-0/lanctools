# The .lanc Format

`lanctools` is built around the `.lanc` file format introduced by Hou et. al.
This is a compressed format which uses a breakpoint-encoded representation
of phased local ancestry data.

Rather than storing ancestry at every variant, each individual is
represented as a sequence of ancestry tracts. Each tract records
the 0-based index where the tract ends (exclusive) and the phased
local ancestry state in that interval.

## Local Ancestry Array

The dense representation corresponds to an array of shape $(n, p, 2)$,
where:

- $n$ is the number of samples
- $p$ is the number of variants
- The final dimension corresponds to the two phased haplotypes (ploidy=2)

The ancestry values in this array are coded as integers corresponding to
the admixing populations in the sample.

## File Format

The first line contains:

`p n`

where $p$ is the number of variants and $n$ is the number of samples.
This is followed by exactly $n$ lines, one per individual.

Each individual line consists of one or more ordered breakpoint records,
separated by spaces. Each breakpoint record has the form:

`<stop>:<anc0><anc1>`

where

- `stop` is the exclusive ending variant index of the current ancestry tract
- `anc0` is the ancestry (integer-coded) at haplotype 0
- `anc1` is the ancestry (integer-coded) at haplotype 1

## Example

The following is a `.lanc` file for 4 variants and 3 samples.

    4 3
    1:00 4:10
    4:11
    1:01 2:11 3:10 4:00

The corresponding dense matrices for haplotype 0 is:

\begin{pmatrix}
0 & 1 & 1 & 1\\
1 & 1 & 1 & 1\\
0 & 1 & 1 & 0
\end{pmatrix}

And for haplotype 1:

\begin{pmatrix}
0 & 0 & 0 & 0\\
1 & 1 & 1 & 1\\
1 & 1 & 0 & 0
\end{pmatrix}

## Limitations

The `.lanc` format has several limitations which may be addressed in future work.

- It can represent a maximum of 10 ancestries
- It cannot represent missing data
- It can only represent autosomes
