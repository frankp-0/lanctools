#include <cstdint>
#include <cstring>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace py = pybind11;

py::tuple query_lanc(
    py::array_t<uint8_t, py::array::c_style | py::array::forcecast> left_haps,
    py::array_t<uint8_t, py::array::c_style | py::array::forcecast> right_haps,
    py::array_t<uint32_t, py::array::c_style | py::array::forcecast>
        breakpoints,
    py::array_t<uint32_t, py::array::c_style | py::array::forcecast> offsets,
    py::array_t<uint32_t, py::array::c_style | py::array::forcecast> indices) {

  const ssize_t n_samples = offsets.shape(0) - 1;
  const ssize_t n_variants = indices.shape(0);

  std::vector<uint8_t> left_buf(left_haps.size());
  std::memcpy(left_buf.data(), left_haps.data(), left_haps.size());

  std::vector<uint8_t> right_buf(right_haps.size());
  std::memcpy(right_buf.data(), right_haps.data(), right_haps.size());

  std::vector<uint32_t> breakpoints_buf(breakpoints.size());
  std::memcpy(breakpoints_buf.data(), breakpoints.data(),
              breakpoints.size() * sizeof(uint32_t));

  std::vector<uint32_t> offsets_buf(offsets.size());
  std::memcpy(offsets_buf.data(), offsets.data(),
              offsets.size() * sizeof(uint32_t));

  std::vector<uint32_t> indices_buf(indices.size());
  std::memcpy(indices_buf.data(), indices.data(),
              indices.size() * sizeof(uint32_t));

  // Allocate outputs
  py::array_t<uint8_t> left_out({n_samples, n_variants});
  py::array_t<uint8_t> right_out({n_samples, n_variants});
  auto left_out_buf = left_out.mutable_unchecked<2>();
  auto right_out_buf = right_out.mutable_unchecked<2>();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (ssize_t i = 0; i < n_samples; ++i) {
    uint32_t start = offsets_buf[i];
    uint32_t end = offsets_buf[i + 1];
    uint32_t end_len = end - start;

    const uint32_t *breakpoints_i = &breakpoints_buf[start];
    const uint8_t *left_i = &left_buf[start];
    const uint8_t *right_i = &right_buf[start];

    uint32_t j = 0;
    for (ssize_t q = 0; q < n_variants; ++q) {
      uint32_t idx = indices_buf[q];

      while (j < end_len && idx >= breakpoints_i[j]) {
        ++j;
      }

      // Clamp to last element
      uint32_t jj = (j == end_len) ? (end_len - 1) : j;

      left_out_buf(i, q) = left_i[jj];
      right_out_buf(i, q) = right_i[jj];
    }
  }

  return py::make_tuple(left_out, right_out);
}

void bind_query_lanc(py::module_ &m) {
  m.def("query_lanc", &query_lanc, "Query local ancestry");
}
