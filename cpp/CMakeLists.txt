# 3.28 is required for FILE_SET headers (3.23), COMPILE_WARNING_AS_ERROR (3.24), SYSTEM (3.25) and EXCLUDE_FROM_ALL
# (3.28) dependencies in FetchContent. Policies are set to NEW up to the maximum version, which is the newest version
# memilio is tested with.
cmake_minimum_required(VERSION 3.28...3.31)

# keep in sync with the latest release on GitHub
project(memilio VERSION 2.4.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# memilio does not use C++20 modules, so skip the dependency scanning that CMake does by default for C++20
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

# Tests, examples and installation default to OFF if memilio is included in another project, e.g. using FetchContent
option(MEMILIO_BUILD_TESTS "Build memilio unit tests." ${PROJECT_IS_TOP_LEVEL})
option(MEMILIO_BUILD_EXAMPLES "Build memilio examples." ${PROJECT_IS_TOP_LEVEL})
option(MEMILIO_BUILD_MODELS "Build memilio models." ON)
option(MEMILIO_BUILD_BENCHMARKS "Build memilio benchmarks with google benchmark." OFF)
option(MEMILIO_ENABLE_INSTALL "Add install rules and the CMake package config for memilio." ${PROJECT_IS_TOP_LEVEL})
option(MEMILIO_USE_BUNDLED_SPDLOG "Download and use spdlog instead of an installed version." ON)
option(MEMILIO_USE_BUNDLED_EIGEN "Download and use Eigen instead of an installed version." ON)
option(MEMILIO_USE_BUNDLED_BOOST "Download and use Boost instead of an installed version." ON)
option(MEMILIO_USE_BUNDLED_JSONCPP "Download and use JsonCpp instead of an installed version." ON)
option(MEMILIO_SANITIZE_ADDRESS "Enable address sanitizer." OFF)
option(MEMILIO_SANITIZE_UNDEFINED "Enable undefined behavior sanitizer." OFF)
option(MEMILIO_BUILD_SHARED_LIBS "Build memilio as a shared library. If OFF, static libraries are built." ON)
option(MEMILIO_ENABLE_MPI "Build memilio with MPI." OFF)
option(MEMILIO_ENABLE_OPENMP "Enable Multithreading with OpenMP." OFF)
option(MEMILIO_ENABLE_HDF5 "Build memilio with HDF5 IO support." ON)
option(MEMILIO_ENABLE_SBML "Build memilio with the SBML importer and imported models." OFF)
option(MEMILIO_ENABLE_WARNINGS "Build memilio with warnings." ON)
option(MEMILIO_ENABLE_WARNINGS_AS_ERRORS "Build memilio with warnings as errors." ON)
option(MEMILIO_ENABLE_IPOPT "Enable numerical optimization with Ipopt, requires a Fortran compiler." OFF)
option(MEMILIO_ENABLE_PROFILING "Enable runtime performance profiling of memilio." OFF)
option(MEMILIO_ENABLE_LIKWID_MARKER "Enable performance measuring with likwid markers." OFF)

cmake_path(CONVERT "${PROJECT_SOURCE_DIR}/.." TO_CMAKE_PATH_LIST PROJECT_ROOT_NORMALIZED NORMALIZE)
set(MEMILIO_BASE_DIR "${PROJECT_ROOT_NORMALIZED}" CACHE PATH
    "Path to the MEmilio project root, used for file I/O. Accessed via mio::base_dir().")
set(MEMILIO_DATA_DIR "${PROJECT_ROOT_NORMALIZED}/data" CACHE PATH
    "Path to the MEmilio data directory. Accessed via mio::data_dir().")

mark_as_advanced(
    MEMILIO_USE_BUNDLED_SPDLOG
    MEMILIO_SANITIZE_ADDRESS
    MEMILIO_SANITIZE_UNDEFINED
    MEMILIO_BUILD_SHARED_LIBS
    MEMILIO_BASE_DIR
    MEMILIO_DATA_DIR
    )

# AppleClang is handled by the Clang options below (all checks use MATCHES "Clang"), but warn about missing support
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message(WARNING "The compiler ID \"AppleClang\" has only limited support. Continuing with \"Clang\" options.")
    add_compile_options(-ffp-model=strict) # to avoid FP differences in CI
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

get_property(MEMILIO_IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(PROJECT_IS_TOP_LEVEL AND NOT MEMILIO_IS_MULTI_CONFIG AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Please choose Release or Debug. Default is Release." FORCE)
endif()

set(BUILD_SHARED_LIBS ${MEMILIO_BUILD_SHARED_LIBS})
if(BUILD_SHARED_LIBS AND MSVC)
    # Current workaround for windows to not build shared libs as we need to plugin export variables
    # TODO: add exports for windows files
    message(WARNING "Building shared libs is not supported on windows. Building static libs instead.")
    set(BUILD_SHARED_LIBS OFF)
endif()

include(GNUInstallDirs)
include(MemilioHelpers)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
# let installed libraries and executables find the installed shared libraries, relative to their own location
# (the build tree RPATH is set separately by CMake)
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path" "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
    set(CMAKE_INSTALL_RPATH "$ORIGIN" "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()

# code coverage analysis
# Note: this only works under linux and with make
# Ninja creates different directory names which do not work together with this scrupt
# as STREQUAL is case-sensitive https://github.com/TriBITSPub/TriBITS/issues/131, also allow DEBUG as accepted input
if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    option(MEMILIO_TEST_COVERAGE "Enable GCov coverage analysis (adds a 'coverage' target)" OFF)
    mark_as_advanced(MEMILIO_TEST_COVERAGE)

    if(MEMILIO_TEST_COVERAGE)
        message(STATUS "Coverage enabled")
        include(CodeCoverage)
        append_coverage_compiler_flags()

        # In addition to standard flags, disable elision and inlining to prevent e.g. closing brackets being marked as
        # uncovered.
        add_compile_options(-fno-elide-constructors -fno-default-inline)
        setup_target_for_coverage_lcov(
            NAME coverage
            EXECUTABLE memilio-test
            LCOV_ARGS --ignore-errors gcov,mismatch
            EXCLUDE "${PROJECT_SOURCE_DIR}/tests*" "${PROJECT_SOURCE_DIR}/examples*"
            "${PROJECT_SOURCE_DIR}/memilio/ad*" "${PROJECT_BINARY_DIR}/*" "/usr*" "${PROJECT_SOURCE_DIR}/*.F"
        )
    endif()
endif()

# set sanitizer compiler flags, only for Debug builds
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    if(MEMILIO_SANITIZE_ADDRESS)
        add_compile_options($<$<CONFIG:Debug>:-fsanitize=address>)
        add_link_options($<$<CONFIG:Debug>:-fsanitize=address>)
    endif()

    if(MEMILIO_SANITIZE_UNDEFINED)
        add_compile_options($<$<CONFIG:Debug>:-fsanitize=undefined>)
        add_link_options($<$<CONFIG:Debug>:-fsanitize=undefined>)
    endif()

    if(MEMILIO_SANITIZE_ADDRESS OR MEMILIO_SANITIZE_UNDEFINED)
        add_compile_options($<$<CONFIG:Debug>:-fno-omit-frame-pointer>)
    endif()
endif()

# define flags to enable most warnings for different compilers
# add flags to each target separately (see memilio_target_warnings) instead of globally so users have the choice to use
# their own flags
set(MEMILIO_CXX_WARNING_FLAGS)
if(MEMILIO_ENABLE_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        list(APPEND MEMILIO_CXX_WARNING_FLAGS /W4)
    else()
        list(APPEND MEMILIO_CXX_WARNING_FLAGS -Wall -Wextra -Wshadow --pedantic)
    endif()
endif()

# exclude some warnings we accept
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    list(APPEND MEMILIO_CXX_WARNING_FLAGS
        -Wno-unknown-warning -Wno-pragmas -Wno-deprecated-copy -Wno-expansion-to-defined)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    list(APPEND MEMILIO_CXX_WARNING_FLAGS
        -Wno-unknown-warning-option -Wno-deprecated -Wno-gnu-zero-variadic-macro-arguments)
endif()

# workarounds for compiler bugs or overzealous optimization/analysis
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    list(APPEND MEMILIO_CXX_WARNING_FLAGS
        $<$<CONFIG:Release>:-Wno-restrict> # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105329
        $<$<CONFIG:Release>:-Wno-maybe-uninitialized> # https://gcc.gnu.org/bugzilla/buglist.cgi?quicksearch=may%20be%20uninitialized
        $<$<CONFIG:Release>:-Wno-stringop-overread>
    )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    list(APPEND MEMILIO_CXX_WARNING_FLAGS $<$<CONFIG:Release>:-Wno-stringop-overread>)
endif()

# warning flags including warnings as errors, for CMake files generated by the SBML importer and memilio-generation
set(MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS ${MEMILIO_CXX_WARNING_FLAGS})
if(MEMILIO_ENABLE_WARNINGS_AS_ERRORS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        list(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS /WX)
    else()
        list(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS -Werror)
    endif()
endif()

# add parts of the project
include(thirdparty/CMakeLists.txt)
add_subdirectory(memilio/ad)
add_subdirectory(memilio)

if(MEMILIO_BUILD_MODELS)
    set(MEMILIO_MODELS
        abm
        d_abm
        ode_secir
        ode_secirts
        ode_secirvvs
        lct_secir
        lct_secir_2_diseases
        glct_secir
        ide_secir
        ide_seir
        ode_seir
        ode_seirv
        ode_seirdb
        ode_seair
        ode_sir
        ode_seir_metapop
        sde_sir
        sde_sirs
        sde_seirvv
        graph_abm
        smm
        hybrid
        ode_mseirs4
    )
    foreach(model IN LISTS MEMILIO_MODELS)
        add_subdirectory(models/${model})
    endforeach()
endif()

if(MEMILIO_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(MEMILIO_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(MEMILIO_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if(MEMILIO_ENABLE_SBML)
    add_subdirectory(sbml_model_generation)
endif()

if(MEMILIO_ENABLE_INSTALL)
    get_property(MEMILIO_MODEL_TARGETS GLOBAL PROPERTY MEMILIO_MODEL_TARGETS)
    install(TARGETS memilio AD ${MEMILIO_MODEL_TARGETS} ${MEMILIO_BUNDLED_TARGETS}
        EXPORT memilio-targets
        FILE_SET HEADERS
    )
    install(EXPORT memilio-targets
        NAMESPACE memilio::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
    )

    include(CMakePackageConfigHelpers)
    configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/memilio-config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/memilio-config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
    )
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/memilio-config-version.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )
    install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/memilio-config-version.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/memilio-config.cmake"
        "${CMAKE_CURRENT_LIST_DIR}/cmake/FindLikwid.cmake"
        DESTINATION
        ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
    )
endif()
