# Fast .lanc Queries

## FlatLanc

Based on the `.lanc` file format, we represent local ancestry in memory using the
[FlatLanc](/docs/api.md#lanctools.FlatLanc) data structure. Rather than storing a
dense local ancestry matrix, `FlatLanc` stores only ancestry transitions (tract
boundaries), which substantially reduces memory usage for typical datasets where
ancestry changes infrequently along the genome.

`FlatLanc` consists of four contiguous one-dimensional arrays:

1. The left haplotype ancestry at each breakpoint, concatenated across all samples.
2. The right haplotype ancestry at each breakpoint, concatenated across all samples.
3. The breakpoint indices marking the end of each ancestry tract, concatenated across all samples.
4. An offset array indicating where each sample's data begins and ends in the three arrays above.

For sample `i`, its data are accessed as

```python
start = offsets[i]
end = offsets[i + 1]

breakpoints_i = breakpoints[start:end]
left_haps_i = left_haps[start:end]
right_haps_i = right_haps[start:end]
```

Each entry of `left_haps_i` and `right_haps_i` corresponds to the ancestry of
the tract ending at the associated breakpoint in `breakpoints_i`.

## Query Algorithm

Given a sorted array of variant indices, ancestry is recovered by simultaneously
scanning the query positions and the sample's ancestry tract.

For each sample:

1. Initialize a pointer to the first ancestry tract.
2. Iterate over the query indices in ascending order.
3. Advance the tract pointer whenever the current query index reaches or passes
   the next breakpoint.
4. Assign the ancestry labels of the current tract to the query position.
5. Continue until all query positions have been processed.

Since both the breakpoint array and the query indices are sorted, the tract
pointer only moves forward and each breakpoint is visited at most once. This is
analogous to the merge step of merge sort or the two-pointer technique commonly
used in sorted-array algorithms.

### Pseudocode

```text
GET_LOCAL_ANCESTRY(
    left_haps,
    right_haps,
    breakpoints,
    offsets,
    indices
)

Require: indices is sorted in ascending order.

Allocate left_out[n_samples][n_variants]
Allocate right_out[n_samples][n_variants]

for each sample do
    start ← offsets[sample]
    end ← offsets[sample + 1]

    end_i ← breakpoints[start:end]
    left_i ← left_haps[start:end]
    right_i ← right_haps[start:end]

    j ← 0

    for each query index idx do
        while j < length(end_i) and idx ≥ end_i[j] do
            j ← j + 1
        end while

        left_out[sample, idx] ← left_i[j]
        right_out[sample, idx] ← right_i[j]
    end for
end for

return left_out, right_out
```

### Complexity

For a sample with **S** ancestry tracts and **V** query positions, the algorithm
runs in **O(S + V)** time because both arrays are traversed only once. Across
**N** samples, the total complexity is

\[
O\left(\sum_{i=1}^{N} S_i + NV\right),
\]

where \(S_i\) is the number of ancestry tracts for sample \(i\).

The outer loop is parallelized with `numba.prange`, allowing independent samples
to be processed concurrently on multicore systems.

An alternative approach is to perform a binary search over the breakpoints for
each query position, requiring $O(V\log S)$ time per sample. For a small number
of query positions this may have a lower constant cost, but for the common case
of querying many variants across many samples, the merge-style scan is more
efficient because each breakpoint is visited at most once per sample and the
sorted query indices are reused across all samples.
