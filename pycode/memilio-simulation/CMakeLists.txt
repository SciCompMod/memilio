cmake_minimum_required(VERSION 3.15)
project(memilio-python LANGUAGES CXX C)

# add in C++ library
set(MEMILIO_CPP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp)
add_subdirectory(${MEMILIO_CPP_DIR} ${CMAKE_CURRENT_BINARY_DIR}/cpp EXCLUDE_FROM_ALL)
set(CMAKE_MODULE_PATH "${MEMILIO_CPP_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module) #help pybind find the right python, fails on multilinux2014 docker images
find_package(pybind11 CONFIG REQUIRED)

# build python extensions
pybind11_add_module(_simulation_secir MODULE
    memilio/simulation/secir.cpp
)
target_link_libraries(_simulation_secir PRIVATE secir)
target_include_directories(_simulation_secir PRIVATE memilio/simulation)
install(TARGETS _simulation_secir LIBRARY DESTINATION memilio)

pybind11_add_module(_simulation MODULE
    memilio/simulation/simulation.cpp
)
target_link_libraries(_simulation PRIVATE memilio)
target_include_directories(_simulation PRIVATE memilio/simulation)
install(TARGETS _simulation LIBRARY DESTINATION memilio)
