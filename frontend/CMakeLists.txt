cmake_minimum_required(VERSION 3.4)

project(epi-frontend VERSION 0.1.0)


if (NOT EMSCRIPTEN)
  message(FATAL_ERROR "This porject must be run with the Emscripten sdk via \"emcmake cmake <cmake_options>\"")
endif(NOT EMSCRIPTEN)

# add in C++ library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cpp ${CMAKE_CURRENT_BINARY_DIR}/cpp EXCLUDE_FROM_ALL)

add_executable(secirjs src/bindings/secir_js.cpp)
target_compile_features(secirjs PRIVATE cxx_std_14)

target_link_libraries(secirjs PRIVATE epidemiology)
set_target_properties(secirjs PROPERTIES LINK_FLAGS  "--bind -s MODULARIZE=1 -s EXPORT_NAME=\"secirjs\"")
