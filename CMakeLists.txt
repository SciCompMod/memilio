cmake_minimum_required(VERSION 3.13)
project(memilio-simulation)

set(CMAKE_CXX_STANDARD "20")
set(CMAKE_CXX_STANDARD_REQUIRED "20")

# add in C++ library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cpp ${CMAKE_CURRENT_BINARY_DIR}/cpp EXCLUDE_FROM_ALL)

add_subdirectory(pycode/memilio-simulation)