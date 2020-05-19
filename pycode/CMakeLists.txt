cmake_minimum_required(VERSION 3.4)
project(epi-py)

option(EPIPY_USE_SYSTEM_PYBIND11 "Use pybind11 installed in the system" OFF)
mark_as_advanced(EPIPY_USE_SYSTEM_PYBIND11)

if(NOT EPIPY_USE_SYSTEM_PYBIND11)
  # Fetch pybind11
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.5
  )
  FetchContent_GetProperties(pybind11)
  if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
  endif()
else()
    find_package(pybind11 REQUIRED)
endif()

# add in C++ library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cpp ${CMAKE_CURRENT_BINARY_DIR}/cpp EXCLUDE_FROM_ALL)

# build secir python extension
pybind11_add_module(_secir MODULE
  epidemiology/secir/secir_py.cpp
)
target_link_libraries(_secir PRIVATE epidemiology)

install(TARGETS _secir DESTINATION epidemiology)
install(FILES epidemiology/secir/secir.py DESTINATION epidemiology)