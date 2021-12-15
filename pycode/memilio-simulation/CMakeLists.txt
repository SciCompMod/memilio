cmake_minimum_required(VERSION 3.10)
project(memilio-python LANGUAGES CXX C)

# add in C++ library
set(MEMILIO_CPP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp)
add_subdirectory(${MEMILIO_CPP_DIR} ${CMAKE_CURRENT_BINARY_DIR}/cpp EXCLUDE_FROM_ALL)
set(CMAKE_MODULE_PATH "${MEMILIO_CPP_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(MEMILIO_USE_SYSTEM_PYBIND11 "Use pybind11 installed on the system instead of Conan package." OFF)

set(_MEMILIO_PYBIND11_VERSION 2.7.1)
if (MEMILIO_USE_CONAN)
    message(STATUS "Memilio: Installing required dependencies for python bindings with Conan.")
    
    include(memilio-conan)

    set(_MEMILIO_REQUIRES "")
    if (NOT MEMILIO_USE_SYSTEM_PYBIND11)
        list(APPEND _MEMILIO_REQUIRES pybind11/${_MEMILIO_PYBIND11_VERSION})
    endif()

    memilio_conan_install(REQUIRES ${_MEMILIO_REQUIRES})

    if (NOT MEMILIO_USE_SYSTEM_PYBIND11)
        set(pybind11_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
    endif()
endif()

find_package(pybind11 CONFIG REQUIRED)

# build python extensions
pybind11_add_module(_simulation_secir MODULE
    memilio/simulation/secir.cpp
)
target_link_libraries(_simulation_secir PRIVATE memilio secir pybind11::pybind11)
target_include_directories(_simulation_secir PRIVATE memilio/simulation)
install(TARGETS _simulation_secir LIBRARY DESTINATION memilio)

pybind11_add_module(_simulation MODULE
    memilio/simulation/simulation.cpp
)
target_link_libraries(_simulation PRIVATE memilio pybind11::pybind11)
target_include_directories(_simulation PRIVATE memilio/simulation)
install(TARGETS _simulation LIBRARY DESTINATION memilio)
