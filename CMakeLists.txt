cmake_minimum_required(VERSION 3.15)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(OpenMP)

pybind11_add_module(_cpp
  src/lanctools/_cpp/bindings.cpp
  src/lanctools/_cpp/flare.cpp
  src/lanctools/_cpp/rfmix.cpp
  src/lanctools/_cpp/lanc.cpp
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(_cpp PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS _cpp DESTINATION ${SKBUILD_PROJECT_NAME})

target_link_libraries(_cpp PRIVATE ZLIB::ZLIB)
