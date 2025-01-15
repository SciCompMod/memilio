cmake_minimum_required(VERSION 3.13)

project(memilio VERSION 1.0.0)

option(MEMILIO_BUILD_TESTS "Build memilio unit tests." ON)
option(MEMILIO_BUILD_EXAMPLES "Build memilio examples." ON)
option(MEMILIO_BUILD_MODELS "Build memilio models." ON)
option(MEMILIO_BUILD_SIMULATIONS "Build memilio simulations that were used for scientific articles." ON)
option(MEMILIO_BUILD_BENCHMARKS "Build memilio benchmarks with google benchmark." OFF)
option(MEMILIO_USE_BUNDLED_SPDLOG "Use spdlog bundled with epi" ON)
option(MEMILIO_USE_BUNDLED_EIGEN "Use eigen bundled with epi" ON)
option(MEMILIO_USE_BUNDLED_BOOST "Use boost bundled with epi (only for epi-io)" ON)
option(MEMILIO_USE_BUNDLED_JSONCPP "Use jsoncpp bundled with epi (only for epi-io)" ON)
option(MEMILIO_SANITIZE_ADDRESS "Enable address sanitizer." OFF)
option(MEMILIO_SANITIZE_UNDEFINED "Enable undefined behavior sanitizer." OFF)
option(MEMILIO_BUILD_SHARED_LIBS "Build memilio as a shared library." ON)
option(MEMILIO_BUILD_STATIC_LIBS "Build memilio as a static library." ON)
option(MEMILIO_ENABLE_MPI "Build memilio with MPI." OFF)
option(MEMILIO_ENABLE_OPENMP "Enable Multithreading with OpenMP." OFF)
option(MEMILIO_ENABLE_WARNINGS "Build memilio with warnings." ON)
option(MEMILIO_ENABLE_WARNINGS_AS_ERRORS "Build memilio with warnings as errors." ON)
option(MEMILIO_ENABLE_IPOPT "Enable numerical optimization with Ipopt, requires a Fortran compiler" OFF)
option(MEMILIO_ENABLE_PROFILING "Enable runtime performance profiling of memilio" OFF)

mark_as_advanced(MEMILIO_USE_BUNDLED_SPDLOG MEMILIO_SANITIZE_ADDRESS MEMILIO_SANITIZE_UNDEFINED)

# try to treat AppleClang as Clang, but warn about missing support 
if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message(WARNING "The compiler ID \"AppleClang\" is not supported, trying to compile with \"Clang\" options.")
    set(CMAKE_CXX_COMPILER_ID "Clang")
endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Please choose Release or Debug. Default is Release." FORCE)
endif()

if(MEMILIO_BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)

    if(MSVC)
        # Current workaround for windows to not build shared libs as we neeed to plugin export variables
        # TODO: add expoerts for windows files
        message(WARNING "Building shared libs is not supported on windows. Building static libs instead.")
        set(BUILD_SHARED_LIBS OFF)
        set(BUILD_STATIC_LIBS ON)
    endif()
else()
    set(BUILD_SHARED_LIBS OFF)
endif()

if(MEMILIO_BUILD_STATIC_LIBS)
    set(BUILD_STATIC_LIBS ON)
else()
    set(BUILD_STATIC_LIBS OFF)
endif()

# Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24 and above
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
endif()

include(GNUInstallDirs) # set to gnu folders. No cache variable so this is not global

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_RPATH "${CMAKE_BINARY_DIR}/lib" "${CMAKE_BINARY_DIR}/bin")

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
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-elide-constructors -fno-default-inline")
        setup_target_for_coverage_lcov(
            NAME coverage
            EXECUTABLE memilio-test
            EXCLUDE "${CMAKE_SOURCE_DIR}/tests*" "${CMAKE_SOURCE_DIR}/simulations*" "${CMAKE_SOURCE_DIR}/examples*"
            "${CMAKE_SOURCE_DIR}/memilio/ad*" "${CMAKE_BINARY_DIR}/*" "/usr*" "${CMAKE_SOURCE_DIR}/*.F"
        )
    endif()
endif()

# set sanitizer compiler flags
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7))
    if(MEMILIO_SANITIZE_ADDRESS)
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address")
        string(APPEND CMAKE_LINKER_FLAGS_DEBUG " -fsanitize=address")
    endif(MEMILIO_SANITIZE_ADDRESS)

    if(MEMILIO_SANITIZE_UNDEFINED)
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=undefined")
        string(APPEND CMAKE_LINKER_FLAGS_DEBUG " -fsanitize=undefined")
    endif(MEMILIO_SANITIZE_UNDEFINED)

    if(MEMILIO_SANITIZE_ADDRESS OR MEMILIO_SANITIZE_UNDEFINED)
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fno-omit-frame-pointer")
        string(APPEND CMAKE_LINKER_FLAGS_DEBUG " -fno-omit-frame-pointer")
    endif(MEMILIO_SANITIZE_ADDRESS OR MEMILIO_SANITIZE_UNDEFINED)
endif((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7))

# define flags to enable most warnings and treat them as errors for different compilers
# add flags to each target separately instead of globally so users have the choice to use their own flags
if(MEMILIO_ENABLE_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
            "/W4;")
    else()
        string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
            "-Wall;-Wextra;-Wshadow;--pedantic;")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
        "-Wno-unknown-warning;-Wno-pragmas;-Wno-deprecated-copy;-Wno-expansion-to-defined;")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
        "-Wno-unknown-warning-option;-Wno-deprecated;-Wno-gnu-zero-variadic-macro-arguments;")
endif()

if(MEMILIO_ENABLE_WARNINGS_AS_ERRORS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
            "/WX")
    else()
        string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
            "-Werror")
    endif()
endif()

# add parts of the project
include(thirdparty/CMakeLists.txt)
add_subdirectory(memilio/ad)
add_subdirectory(memilio)

if(MEMILIO_BUILD_MODELS)
    add_subdirectory(models/abm)
    add_subdirectory(models/ode_secir)
    add_subdirectory(models/ode_secirts)
    add_subdirectory(models/ode_secirvvs)
    add_subdirectory(models/lct_secir)
    add_subdirectory(models/glct_secir)
    add_subdirectory(models/ide_secir)
    add_subdirectory(models/ide_seir)
    add_subdirectory(models/ode_seir)
    add_subdirectory(models/ode_seair)
    add_subdirectory(models/ode_sir)
    add_subdirectory(models/sde_sir)
    add_subdirectory(models/sde_sirs)
    add_subdirectory(models/sde_seirvv)
    add_subdirectory(models/graph_abm)
endif()

if(MEMILIO_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(MEMILIO_BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(MEMILIO_BUILD_SIMULATIONS)
    add_subdirectory(simulations)
endif()

if(MEMILIO_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

install(TARGETS memilio
    EXPORT memilio-targets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY memilio DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN memilio/*/*.h)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/memilio DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN memilio/*/*.h)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/memilio-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/memilio-config.cmake
    INSTALL_DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/memilio-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
    "${CMAKE_CURRENT_BINARY_DIR}/memilio-config-version.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/memilio-config.cmake"
    DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
)
