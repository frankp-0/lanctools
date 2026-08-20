#pragma once
#include <pybind11/pybind11.h>

namespace py = pybind11;

/// Parse an RFMix MSP file into column-oriented ancestry tract data.
///
/// The returned dictionary contains the columns `sample`, `chrom`, `spos`,
/// `epos`, `anc0`, and `anc1`, which are consumed by the Python conversion
/// pipeline.
py::dict read_rfmix(const std::string &msp_file);

/// Add the RFMix reader to the Python extension module.
void bind_rfmix(py::module_ &m);
