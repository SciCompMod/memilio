cmake_minimum_required(VERSION 3.4)

project(epi-frontend VERSION 0.1.0)

if (NOT EMSCRIPTEN)
  message(FATAL_ERROR "This project must be run with the Emscripten sdk via \"emcmake cmake <cmake_options>\"")
endif(NOT EMSCRIPTEN)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/public)

# add in C++ library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../cpp ${CMAKE_CURRENT_BINARY_DIR}/cpp EXCLUDE_FROM_ALL)

add_executable(secirjs src/bindings/secir_js.cpp)
target_compile_features(secirjs PRIVATE cxx_std_14)

target_link_libraries(secirjs PRIVATE epidemiology)
set_target_properties(secirjs PROPERTIES LINK_FLAGS  "--bind -s MODULARIZE=1 -s EXPORT_NAME=\"secirjs\" -s INITIAL_MEMORY=536870912")

if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s ASSERTIONS=1 -O0 -g3 -fsanitize=address")
    add_definitions(-DDEBUG)
    add_link_options("SHELL:-s DISABLE_EXCEPTION_CATCHING=0")
    add_link_options("SHELL:-s ERROR_ON_UNDEFINED_SYMBOLS=1")
    add_link_options("SHELL:-s DEMANGLE_SUPPORT=1")
endif()

install(FILES ${CMAKE_BINARY_DIR}/secirjs.js ${CMAKE_BINARY_DIR}/secirjs.wasm DESTINATION .)