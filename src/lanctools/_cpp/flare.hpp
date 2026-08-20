#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

/// Parse a gzipped FLARE ancestry VCF into column-oriented tract data.
///
/// The returned dictionary contains the columns `sample`, `chrom`, `spos`,
/// `epos`, `anc0`, and `anc1`. Coordinates are genomic positions; the Python
/// layer maps them to PLINK variant indices before writing a `.lanc` file.
py::dict read_flare(const std::string &flare_file);

/// Add the FLARE reader to the Python extension module.
void bind_flare(py::module_ &m);
