#include "flare.hpp"
#include "lanc.hpp"
#include "rfmix.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(_cpp, m) {
  m.doc() = "lanctools C++ backend";
  bind_rfmix(m);
  bind_flare(m);
  bind_query_lanc(m);
}
