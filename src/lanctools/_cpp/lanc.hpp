#pragma once

#include <cstdint>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

py::tuple query_lanc(
    py::array_t<uint8_t, py::array::c_style | py::array::forcecast> left_haps,
    py::array_t<uint8_t, py::array::c_style | py::array::forcecast> right_haps,
    py::array_t<uint32_t, py::array::c_style | py::array::forcecast>
        breakpoints,
    py::array_t<uint32_t, py::array::c_style | py::array::forcecast> offsets,
    py::array_t<uint32_t, py::array::c_style | py::array::forcecast> indices);
void bind_query_lanc(py::module_ &m);
