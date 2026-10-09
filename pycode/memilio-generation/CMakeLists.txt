cmake_minimum_required(VERSION 3.28)
project(memilio-python)

include(GNUInstallDirs) # set to gnu folders. No cache variable so this is not global

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# add in C++ library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp ${CMAKE_CURRENT_BINARY_DIR}/cpp EXCLUDE_FROM_ALL)

# libclang is usually installed next to the clang compiler
cmake_path(GET CMAKE_CXX_COMPILER PARENT_PATH COMPILER_BIN_DIR)
find_library(LIBCLANG_PATH clang HINTS "${COMPILER_BIN_DIR}/../lib")

# config.json is read by the generation scripts from the package directory, so it has to be generated in the source tree.
# It contains local paths and is ignored by git.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/memilio/tools/config.json.txt ${CMAKE_CURRENT_SOURCE_DIR}/memilio/tools/config.json)
