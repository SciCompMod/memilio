cmake_minimum_required(VERSION 3.13)

project(memilio VERSION 1.0.0)

set(CMAKE_CXX_STANDARD "20")
set(CMAKE_CXX_STANDARD_REQUIRED "20")

option(MEMILIO_BUILD_TESTS "Build memilio unit tests." ON)
option(MEMILIO_BUILD_EXAMPLES "Build memilio examples." ON)
option(MEMILIO_BUILD_MODELS "Build memilio models." ON)
option(MEMILIO_BUILD_BENCHMARKS "Build memilio benchmarks with google benchmark." OFF)
option(MEMILIO_BUILD_SBML_MODELS "Build the SBML importer and imported models." ON)
option(MEMILIO_USE_BUNDLED_SPDLOG "Use spdlog bundled with epi" ON)
option(MEMILIO_USE_BUNDLED_EIGEN "Use eigen bundled with epi" ON)
option(MEMILIO_USE_BUNDLED_BOOST "Use boost bundled with epi (only for epi-io)" ON)
option(MEMILIO_USE_BUNDLED_JSONCPP "Use jsoncpp bundled with epi (only for epi-io)" ON)
option(MEMILIO_SANITIZE_ADDRESS "Enable address sanitizer." OFF)
option(MEMILIO_SANITIZE_UNDEFINED "Enable undefined behavior sanitizer." OFF)
option(MEMILIO_BUILD_SHARED_LIBS "Build memilio as a shared library." ON)
option(MEMILIO_BUILD_STATIC_LIBS "Build memilio as a static library." ON)
option(MEMILIO_ENABLE_MPI "Build memilio with MPI." ON)
option(MEMILIO_ENABLE_OPENMP "Enable Multithreading with OpenMP." OFF)
option(MEMILIO_ENABLE_WARNINGS "Build memilio with warnings." ON)
option(MEMILIO_ENABLE_WARNINGS_AS_ERRORS "Build memilio with warnings as errors." ON)
option(MEMILIO_ENABLE_IPOPT "Enable numerical optimization with Ipopt, requires a Fortran compiler." OFF)
option(MEMILIO_ENABLE_PROFILING "Enable runtime performance profiling of memilio." OFF)
option(MEMILIO_ENABLE_LIKWID_MARKER "Enable performance measuring with likwid markers." OFF)

mark_as_advanced(MEMILIO_USE_BUNDLED_SPDLOG MEMILIO_SANITIZE_ADDRESS MEMILIO_SANITIZE_UNDEFINED)

# try to treat AppleClang as Clang, but warn about missing support
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
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

file(TO_CMAKE_PATH "${PROJECT_SOURCE_DIR}/.." MEMILIO_BASE_DIR)

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
            LCOV_ARGS --ignore-errors gcov,mismatch
            EXCLUDE "${CMAKE_SOURCE_DIR}/tests*" "${CMAKE_SOURCE_DIR}/examples*"
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

# exclude some warnings we accept
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
        "-Wno-unknown-warning;-Wno-pragmas;-Wno-deprecated-copy;-Wno-expansion-to-defined;")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
        "-Wno-unknown-warning-option;-Wno-deprecated;-Wno-gnu-zero-variadic-macro-arguments;")
endif()

# woyrkarounds for compiler bugs or overzealous optimization/analysis
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RELEASE")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
            "-Wno-restrict;" # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105329
            "-Wno-maybe-uninitialized;" # https://gcc.gnu.org/bugzilla/buglist.cgi?quicksearch=may%20be%20uninitialized
            "-Wno-stringop-overread;"
        )
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        string(APPEND MEMILIO_CXX_FLAGS_ENABLE_WARNING_ERRORS
            "-Wno-stringop-overread;"
        )
    endif()
endif()

# finalize string by setting warnings as errors
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
    add_subdirectory(models/d_abm)
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
    add_subdirectory(models/smm)
    add_subdirectory(models/hybrid)
    add_subdirectory(models/ode_mseirs4)
endif()

if(MEMILIO_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(MEMILIO_BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(MEMILIO_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if(MEMILIO_BUILD_SBML_MODELS)
    add_subdirectory(sbml_model_generation)
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
