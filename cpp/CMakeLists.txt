cmake_minimum_required(VERSION 3.24)

project(memilio 
    VERSION 2.0.0 
    DESCRIPTION "A library for the simulations of infectious dynamics."
    HOMEPAGE_URL "https://github.com/SciCompMod/memilio" 
    LANGUAGES CXX)

#-------------------------------------------------------
# System requirements and compatibility
#-------------------------------------------------------
# Define minimum supported compiler versions
set(MEMILIO_MIN_GCC_VERSION "9.0")
set(MEMILIO_MIN_CLANG_VERSION "6.0") 
set(MEMILIO_MIN_MSVC_VERSION "19.20")

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type: Debug, Release, RelWithDebInfo, MinSizeRel" FORCE)
    message(STATUS "No build type specified. Setting to Release.")
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 

#-------------------------------------------------------
# Build options
#-------------------------------------------------------
# Library options
option(MEMILIO_ENABLE_INSTALL "Enable installation of memilio" OFF)
option(MEMILIO_BUILD_SHARED_LIBS "Build memilio as a shared library (faster compilation, useful for development)" OFF)
if(MEMILIO_BUILD_SHARED_LIBS AND MSVC)
    message(WARNING "Building shared libs is not supported on Windows. Building static libs instead.")
    set(MEMILIO_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
endif()
set(BUILD_SHARED_LIBS ${MEMILIO_BUILD_SHARED_LIBS})

# Component options
option(MEMILIO_BUILD_TESTS "Build memilio unit tests." ON)
option(MEMILIO_BUILD_EXAMPLES "Build memilio examples." ON)
option(MEMILIO_BUILD_MODELS "Build memilio models." ON)
option(MEMILIO_BUILD_BENCHMARKS "Build memilio benchmarks with google benchmark." ON)

# Third-party dependency options
option(MEMILIO_USE_BUNDLED_EIGEN "Use eigen bundled with memilio" ON)
option(MEMILIO_USE_BUNDLED_BOOST "Use boost bundled with memilio (only for memilio-io)" ON)
option(MEMILIO_USE_BUNDLED_JSONCPP "Use jsoncpp bundled with memilio (only for memilio-io)" ON)
option(MEMILIO_USE_BUNDLED_SPDLOG "Use spdlog bundled with memilio" ON)

# Performance and debugging options
option(MEMILIO_ENABLE_MPI "Build memilio with MPI." OFF)
option(MEMILIO_ENABLE_OPENMP "Enable Multithreading with OpenMP." OFF)
option(MEMILIO_ENABLE_IPOPT "Enable numerical optimization with Ipopt, requires a Fortran compiler" OFF)
option(MEMILIO_ENABLE_PROFILING "Enable runtime performance profiling of memilio" OFF) 
option(MEMILIO_SANITIZE_ADDRESS "Enable address sanitizer (Debug builds only)." OFF)
option(MEMILIO_SANITIZE_UNDEFINED "Enable undefined behavior sanitizer (Debug builds only)." OFF)
option(MEMILIO_TEST_COVERAGE "Enable GCov coverage analysis (adds a 'coverage' target)" OFF)

# Build quality options
option(MEMILIO_ENABLE_WARNINGS "Build memilio with warnings." ON)
option(MEMILIO_ENABLE_WARNINGS_AS_ERRORS "Build memilio with warnings as errors." ON)

mark_as_advanced(
    # Hide advanced build options
    MEMILIO_SANITIZE_ADDRESS
    MEMILIO_SANITIZE_UNDEFINED
    MEMILIO_TEST_COVERAGE
    MEMILIO_ENABLE_PROFILING
    MEMILIO_USE_BUNDLED_EIGEN
    MEMILIO_USE_BUNDLED_BOOST
    MEMILIO_USE_BUNDLED_JSONCPP
    MEMILIO_USE_BUNDLED_SPDLOG
    MEMILIO_ENABLE_IPOPT
    MEMILIO_ENABLE_INSTALL
)
file(TO_CMAKE_PATH "${PROJECT_SOURCE_DIR}/.." MEMILIO_BASE_DIR)
#-------------------------------------------------------
# Helper functions
#-------------------------------------------------------

# Function to check compiler version against minimum
function(check_compiler_version COMPILER_ID CURRENT_VERSION MIN_VERSION)
    if(CMAKE_CXX_COMPILER_ID STREQUAL ${COMPILER_ID} AND 
       CURRENT_VERSION VERSION_LESS ${MIN_VERSION})
        message(WARNING "${COMPILER_ID} version ${CURRENT_VERSION} is below the recommended version ${MIN_VERSION}. Some features may not work correctly.")
    endif()
endfunction()

# Function to add a component conditionally
function(add_memilio_component COMPONENT_TYPE)
    string(TOUPPER ${COMPONENT_TYPE} COMPONENT_TYPE_UPPER)
    if(MEMILIO_BUILD_${COMPONENT_TYPE_UPPER})
        message(STATUS "Building memilio ${COMPONENT_TYPE}")
        add_subdirectory(${COMPONENT_TYPE})
    endif()
endfunction()

# Function to install header directories
function(install_headers DIR)
    install(
        DIRECTORY ${DIR}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING
            PATTERN "*.h"
            PATTERN "*.hpp"
    )
endfunction()

# Function to print a configuration item with consistent formatting
function(print_config_item ITEM_NAME ITEM_VALUE)
    string(LENGTH "${ITEM_NAME}" NAME_LEN)
    set(PADDING_LEN 25)
    math(EXPR SPACES_NEEDED "${PADDING_LEN} - ${NAME_LEN}")    
    string(REPEAT " " ${SPACES_NEEDED} PADDING_SPACES)
    if(SPACES_NEEDED LESS 0)
        set(PADDING_SPACES "")
    endif()    
    message(STATUS "  ${ITEM_NAME}:${PADDING_SPACES} ${ITEM_VALUE}")
endfunction()

#-------------------------------------------------------
# Check for required dependencies
#-------------------------------------------------------
if(MEMILIO_ENABLE_MPI)
    find_package(MPI QUIET)
    if(NOT MPI_FOUND)
        message(FATAL_ERROR "MPI not found. MEMILIO_ENABLE_MPI will be disabled.")
        set(MEMILIO_ENABLE_MPI OFF CACHE BOOL "Disabled because MPI was not found" FORCE)
        return()
    endif()
    message(STATUS "MPI found: ${MPI_VERSION}")
endif()

message(STATUS "MEMILIO_ENABLE_IPOPT is: ${MEMILIO_ENABLE_IPOPT}")
if(MEMILIO_ENABLE_IPOPT)
    enable_language(Fortran)
    message(STATUS "Ipopt enabled, ensuring Fortran language support is enabled.")
endif()

if(MEMILIO_ENABLE_OPENMP)
    find_package(OpenMP QUIET) # Find quietly first
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found via find_package.")
        set(MEMILIO_OPENMP_FOUND TRUE)
    elseif(APPLE)
        message(STATUS "OpenMP not found by find_package, attempting macOS Homebrew workaround...")
        # (Homebrew lookup logic remains the same)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64" OR APPLE_SILICON)
            set(HOMEBREW_PREFIX "/opt/homebrew")
        else()
            set(HOMEBREW_PREFIX "/usr/local")
        endif()
        if(EXISTS "${HOMEBREW_PREFIX}/opt/libomp/lib/cmake/libomp")
             list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}/opt/libomp")
             find_package(OpenMP QUIET) # Try finding again with the hint
             if (OpenMP_CXX_FOUND)
                 message(STATUS "Found OpenMP via Homebrew libomp CMake package.")
                 set(MEMILIO_OPENMP_FOUND TRUE)
             endif()
        endif()
        # Manual target creation fallback if CMake package didn't work/exist
        if (NOT MEMILIO_OPENMP_FOUND AND
            EXISTS "${HOMEBREW_PREFIX}/opt/libomp/include" AND
            EXISTS "${HOMEBREW_PREFIX}/opt/libomp/lib/libomp.dylib")
                 message(STATUS "Found Homebrew libomp files, attempting manual target setup for OpenMP.")
                 # Ensure target doesn't exist before creating
                 if(NOT TARGET OpenMP::OpenMP_CXX)
                     add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
                     set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                         INTERFACE_COMPILE_OPTIONS "-Xpreprocessor -fopenmp" # Clang flag
                         INTERFACE_INCLUDE_DIRECTORIES "${HOMEBREW_PREFIX}/opt/libomp/include"
                         INTERFACE_LINK_LIBRARIES "${HOMEBREW_PREFIX}/opt/libomp/lib/libomp.dylib")
                 endif()
                 set(MEMILIO_OPENMP_FOUND TRUE) # Mark as found via manual setup
        endif()
    endif()

    # Final check and disable if not found
    if(NOT MEMILIO_OPENMP_FOUND)
        message(WARNING "OpenMP not found or configured. Install OpenMP (e.g., 'brew install libomp' on macOS) or disable MEMILIO_ENABLE_OPENMP.")
        set(MEMILIO_ENABLE_OPENMP OFF CACHE BOOL "Disabled because OpenMP was not found/configured" FORCE)
    endif()
endif()

#-------------------------------------------------------
# Compiler detection and platform-specific settings
#-------------------------------------------------------
# Detect compiler and set appropriate flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(STATUS "Detected MSVC compiler (${CMAKE_CXX_COMPILER_VERSION})")
    set(USING_MSVC TRUE)
    check_compiler_version("MSVC" ${CMAKE_CXX_COMPILER_VERSION} ${MEMILIO_MIN_MSVC_VERSION})
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Detected GNU compiler (${CMAKE_CXX_COMPILER_VERSION})")
    set(USING_GNU TRUE)
    check_compiler_version("GNU" ${CMAKE_CXX_COMPILER_VERSION} ${MEMILIO_MIN_GCC_VERSION})
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message(STATUS "Detected AppleClang compiler (${CMAKE_CXX_COMPILER_VERSION} - not officially supported")
    set(USING_CLANG_COMPATIBLE TRUE)
    set(USING_APPLECLANG TRUE)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
        message(STATUS "Detected Apple Silicon (ARM64) architecture")
        set(APPLE_SILICON TRUE)
        set(APPLE_SILICON_FLAGS "-march=native")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "Detected Clang compiler (${CMAKE_CXX_COMPILER_VERSION})")
    set(USING_CLANG_COMPATIBLE TRUE)
    check_compiler_version("Clang" ${CMAKE_CXX_COMPILER_VERSION} ${MEMILIO_MIN_CLANG_VERSION}) 
else()
    message(WARNING "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}. Build may not work correctly.")
endif()

# Set platform-specific settings
if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif()
if(APPLE_SILICON)
    add_compile_options(${APPLE_SILICON_FLAGS})
    message(STATUS "Applied Apple Silicon optimization flags: ${APPLE_SILICON_FLAGS}")
endif()

#-------------------------------------------------------
# Basic build system configuration
#-------------------------------------------------------
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(GNUInstallDirs) # Provides standard installation directories

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")

# Generate version header
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/memilio/version.h
)

#-------------------------------------------------------
# Code coverage configuration
#-------------------------------------------------------
set(MEMILIO_COVERAGE_FLAGS "")
if((CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG") AND MEMILIO_TEST_COVERAGE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "Coverage enabled (GCC/Clang Debug build)")
        include(CodeCoverage) 

        # Define flags for coverage generation and better reporting
        set(MEMILIO_COVERAGE_FLAGS "--coverage;-fno-elide-constructors;-fno-default-inline;-fno-inline") 

        # Define standard exclusions for coverage reports (used by setup_target_for_coverage_lcov)
        set(COVERAGE_EXCLUDES
            "${CMAKE_SOURCE_DIR}/tests*"
            "${CMAKE_SOURCE_DIR}/simulations*"
            "${CMAKE_SOURCE_DIR}/examples*"
            "${CMAKE_SOURCE_DIR}/memilio/ad*"
            "${CMAKE_BINARY_DIR}/*"
            "/usr*"
            "${CMAKE_SOURCE_DIR}/*.F" # Assuming Fortran files if any
            "${CMAKE_SOURCE_DIR}/thirdparty/*" # Exclude fetched sources
            "${PROJECT_BINARY_DIR}/*" # Exclude build artifacts more broadly
        )

        set(MEMILIO_COVERAGE_SETUP_NEEDED TRUE)

        if(CMAKE_GENERATOR MATCHES "Ninja")
            message(WARNING "Using Ninja generator with code coverage may affect results.")
        endif()
    else()
        message(WARNING "Coverage requested but not supported for compiler ${CMAKE_CXX_COMPILER_ID} in this script. Disabling.")
        set(MEMILIO_TEST_COVERAGE OFF CACHE BOOL "Disabled: Compiler not supported" FORCE)
    endif()
endif()

# Function to apply coverage flags (call this in tests/CMakeLists.txt)
function(memilio_set_coverage_flags target_name)
    if(MEMILIO_TEST_COVERAGE AND MEMILIO_COVERAGE_FLAGS)
        # Apply flags using generator expression for Debug config only
        target_compile_options(${target_name} PRIVATE $<$<CONFIG:Debug>:${MEMILIO_COVERAGE_FLAGS}>)
        target_link_options(${target_name} PRIVATE $<$<CONFIG:Debug>:--coverage>) # Linker flag often needed too
    endif()
endfunction()

#-------------------------------------------------------
# Sanitizer configuration
#-------------------------------------------------------
set(MEMILIO_SANITIZER_COMPILE_FLAGS "")
set(MEMILIO_SANITIZER_LINK_FLAGS "")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" OR USING_CLANG_COMPATIBLE) # Check using refined variables
    set(SANITIZERS_ENABLED FALSE)
    if(MEMILIO_SANITIZE_ADDRESS)
        list(APPEND MEMILIO_SANITIZER_COMPILE_FLAGS "-fsanitize=address")
        list(APPEND MEMILIO_SANITIZER_LINK_FLAGS "-fsanitize=address")
        set(SANITIZERS_ENABLED TRUE)
    endif()
    if(MEMILIO_SANITIZE_UNDEFINED)
        list(APPEND MEMILIO_SANITIZER_COMPILE_FLAGS "-fsanitize=undefined")
        list(APPEND MEMILIO_SANITIZER_LINK_FLAGS "-fsanitize=undefined")
        set(SANITIZERS_ENABLED TRUE)
    endif()

    if(SANITIZERS_ENABLED)
        list(APPEND MEMILIO_SANITIZER_COMPILE_FLAGS "-fno-omit-frame-pointer")
        # Linker flags already have the sanitizer flags
        message(STATUS "Sanitizers enabled for Debug builds: ${MEMILIO_SANITIZER_COMPILE_FLAGS}")
    endif()
else()
    if(MEMILIO_SANITIZE_ADDRESS OR MEMILIO_SANITIZE_UNDEFINED)
        message(WARNING "Sanitizers requested but not supported for compiler ${CMAKE_CXX_COMPILER_ID}.")
    endif()
endif()

# Function to apply sanitizer flags (call this after defining targets)
function(memilio_set_sanitizer_flags target_name)
    if(MEMILIO_SANITIZER_COMPILE_FLAGS)
        # Apply flags using generator expression for Debug config only
        target_compile_options(${target_name} PRIVATE $<$<CONFIG:Debug>:${MEMILIO_SANITIZER_COMPILE_FLAGS}>)
        target_link_options(${target_name} PRIVATE $<$<CONFIG:Debug>:${MEMILIO_SANITIZER_LINK_FLAGS}>)
    endif()
endfunction()




#-------------------------------------------------------
# Warning and compiler flags configuration
#-------------------------------------------------------
# Define flags to enable most warnings and treat them as errors for different compilers.
# Initialize variables for different warning configurations
set(MEMILIO_CXX_FLAGS_ENABLE_WARNINGS "")
set(MEMILIO_CXX_FLAGS_DISABLE_WARNINGS "")
set(MEMILIO_CXX_FLAGS_WARNINGS_AS_ERRORS "")

# Set compiler-specific warning flags
if(MEMILIO_ENABLE_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        # MSVC warning flags
        set(MEMILIO_CXX_FLAGS_ENABLE_WARNINGS "/W4;")
    else()
        # GCC/Clang common warning flags
        set(MEMILIO_CXX_FLAGS_ENABLE_WARNINGS 
            "-Wall" 
            "-Wextra" 
            "-Wshadow" 
            "--pedantic"
        )
        
        # Compiler-specific warning suppressions
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            # GCC-specific suppressions
            list(APPEND MEMILIO_CXX_FLAGS_DISABLE_WARNINGS
                "-Wno-unknown-warning"
                "-Wno-pragmas"
                "-Wno-deprecated-copy"
                "-Wno-expansion-to-defined"
            )
        elseif(USING_CLANG_COMPATIBLE)
            # Clang/AppleClang suppressions
            list(APPEND MEMILIO_CXX_FLAGS_DISABLE_WARNINGS
                "-Wno-unknown-warning-option"
                "-Wno-deprecated"
                "-Wno-gnu-zero-variadic-macro-arguments"
            )
        endif()
    endif()

    # Set warnings-as-errors flags if enabled
    if(MEMILIO_ENABLE_WARNINGS_AS_ERRORS)
        if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            set(MEMILIO_CXX_FLAGS_WARNINGS_AS_ERRORS "/WX")
        else()
            set(MEMILIO_CXX_FLAGS_WARNINGS_AS_ERRORS "-Werror")
        endif()
    endif()
endif()

# Combine all warning flags into one list
set(MEMILIO_CXX_WARNING_FLAGS
    ${MEMILIO_CXX_FLAGS_ENABLE_WARNINGS}
    ${MEMILIO_CXX_FLAGS_DISABLE_WARNINGS}
    ${MEMILIO_CXX_FLAGS_WARNINGS_AS_ERRORS}
)

# Define a function to apply warning flags to targets
function(memilio_set_warning_flags target_name)
    if(MEMILIO_ENABLE_WARNINGS)
        target_compile_options(${target_name} PRIVATE ${MEMILIO_CXX_WARNING_FLAGS})
    endif()
endfunction()



#-------------------------------------------------------
# Add project components
#-------------------------------------------------------
# First add third-party dependencies
include(FetchContent)
include(thirdparty/CMakeLists.txt)

# Core library 
add_subdirectory(memilio)

# Conditional components based on build options
if(MEMILIO_BUILD_MODELS)
    message(STATUS "Building memilio models")
    add_subdirectory(models)
endif()

# Use the component function for other directories
add_memilio_component(examples)
add_memilio_component(tests)
add_memilio_component(simulations)
add_memilio_component(benchmarks)


#-------------------------------------------------------
# Package configuration
#-------------------------------------------------------
include(CMakePackageConfigHelpers)

# Create package config files
configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/memilio-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/memilio-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
    PATH_VARS 
    CMAKE_INSTALL_LIBDIR
    CMAKE_INSTALL_INCLUDEDIR
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/memilio-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
    )

#-------------------------------------------------------
# Installation configuration
#-------------------------------------------------------
# Check if installation is enabled
if(MEMILIO_ENABLE_INSTALL)
    message(STATUS "Installation enabled.")

    # Set installation prefix if not set
    if(NOT CMAKE_INSTALL_PREFIX)
        set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/install" CACHE PATH "Installation directory" FORCE)
        message(STATUS "No installation prefix specified. Setting to ${CMAKE_INSTALL_PREFIX}.")
    endif()

    # Set default install component to runtime
    set(CMAKE_INSTALL_COMPONENT_DEFAULT "runtime")

    install(TARGETS memilio
        # EXPORT memilio-targets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT development
    )

    # Install header files using the helper function
    install_headers(memilio)
    install_headers(${CMAKE_CURRENT_BINARY_DIR}/include)

    # Export targets for use in client projects
    install(
        EXPORT memilio-targets
        FILE memilio-targets.cmake
        NAMESPACE memilio::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
        COMPONENT Development
    )

    # Install CMake config files
    install(
        FILES
            "${CMAKE_CURRENT_BINARY_DIR}/memilio-config-version.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/memilio-config.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/memilio
    )
endif()


#-------------------------------------------------------
# Configuration summary
#-------------------------------------------------------
message(STATUS "")
message(STATUS "-------------------- Memilio Configuration Summary --------------------")
print_config_item("Version" "${PROJECT_VERSION}")
print_config_item("Build type" "${CMAKE_BUILD_TYPE}")
print_config_item("Generator" "${CMAKE_GENERATOR}")
print_config_item("Installaion" "${MEMILIO_ENABLE_INSTALL} (Prefix: ${CMAKE_INSTALL_PREFIX})")

message(STATUS "")
message(STATUS "-- Compiler Info --")
print_config_item("C++ Compiler ID" "${CMAKE_CXX_COMPILER_ID}")
print_config_item("C++ Compiler" "${CMAKE_CXX_COMPILER}")
print_config_item("C++ Version" "${CMAKE_CXX_COMPILER_VERSION}")
print_config_item("C++ Standard" "C++${CMAKE_CXX_STANDARD}")
if(MEMILIO_ENABLE_IPOPT AND CMAKE_Fortran_COMPILER)
    print_config_item("Fortran Compiler ID" "${CMAKE_Fortran_COMPILER_ID}")
    print_config_item("Fortran Compiler" "${CMAKE_Fortran_COMPILER}")
    print_config_item("Fortran Version" "${CMAKE_Fortran_COMPILER_VERSION}")
endif()

message(STATUS "")
message(STATUS "-- Components --")
print_config_item("Core Library" "ON") # Assuming memilio target exists
print_config_item("Models" "${MEMILIO_BUILD_MODELS}")
print_config_item("Examples" "${MEMILIO_BUILD_EXAMPLES}")
print_config_item("Tests" "${MEMILIO_BUILD_TESTS}")
print_config_item("Benchmarks" "${MEMILIO_BUILD_BENCHMARKS}")

message(STATUS "")
message(STATUS "-- Dependencies & Features --")
# Check third-party library status using targets or MEMILIO_HAS_* variables
print_config_item("Eigen3" "${MEMILIO_HAS_EIGEN}")
print_config_item("spdlog" "${MEMILIO_HAS_SPDLOG}")
print_config_item("Boost" "${MEMILIO_HAS_EIGEN} ")
print_config_item("HDF5" "${MEMILIO_HAS_HDF5}")
print_config_item("JsonCpp" "${MEMILIO_HAS_JSONCPP}")
print_config_item("Ipopt" "${IPOPT_FOUND} (Enabled: ${MEMILIO_ENABLE_IPOPT})")
print_config_item("gperftools" "${GPERFTOOLS_FOUND} (Enabled: ${MEMILIO_ENABLE_PROFILING})")
print_config_item("OpenMP" "${MEMILIO_OPENMP_FOUND} (Enabled: ${MEMILIO_ENABLE_OPENMP})")
print_config_item("MPI" "${MPI_FOUND} (Enabled: ${MEMILIO_ENABLE_MPI})")

message(STATUS "")
message(STATUS "-- Build Flags --")
print_config_item("Shared libs" "${MEMILIO_BUILD_SHARED_LIBS}")
print_config_item("Warnings" "${MEMILIO_ENABLE_WARNINGS} (As Errors: ${MEMILIO_ENABLE_WARNINGS_AS_ERRORS})")
set(SANITIZER_SUMMARY "OFF")
if(MEMILIO_SANITIZER_COMPILE_FLAGS) 
    set(SANITIZER_SUMMARY "ON (Debug Only)") 
endif()
print_config_item("Sanitizers (ASan/UBSan)" "${SANITIZER_SUMMARY}")
set(COVERAGE_SUMMARY "OFF")
if(MEMILIO_TEST_COVERAGE) 
    set(COVERAGE_SUMMARY "ON (Debug Only, GCC/Clang)") 
endif()
print_config_item("Coverage (gcov)" "${COVERAGE_SUMMARY}")
message(STATUS "======================================================================")