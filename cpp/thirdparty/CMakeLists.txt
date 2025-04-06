include(FetchContent)

#-------------------------------------------------------
# Bundled Library Versions
#-------------------------------------------------------
set(MEMILIO_BUNDLED_EIGEN_VERSION "3.4.0")
set(MEMILIO_BUNDLED_SPDLOG_VERSION "1.15.0")
set(MEMILIO_BUNDLED_BOOST_VERSION "1.84.0")
set(MEMILIO_MINIMAL_BOOST_VERSION "1.76.0") # For system find_package
set(MEMILIO_BUNDLED_JSONCPP_VERSION "1.9.6")
set(MEMILIO_BUNDLED_RANDOM123_VERSION "v1.14.0")
set(MEMILIO_BUNDLED_IPOPT_VERSION "3.14.12")

#-------------------------------------------------------
# Gperftools (Profiling - Optional, Unix-only)
#-------------------------------------------------------
if(MEMILIO_ENABLE_PROFILING)
    if(NOT UNIX)
        message(FATAL_ERROR "Profiling requires a Unix-like system.")
    endif()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GPERFTOOLS REQUIRED libprofiler) # Find using pkg-config

    if(GPERFTOOLS_FOUND)
        # Create an imported target for gperftools
        add_library(memilio::gperftools INTERFACE IMPORTED)
        target_include_directories(memilio::gperftools SYSTEM INTERFACE ${GPERFTOOLS_INCLUDE_DIRS})
        target_link_libraries(memilio::gperftools INTERFACE ${GPERFTOOLS_LINK_LIBRARIES})
        message(STATUS "Profiling enabled: Found system gperftools.")
    else()
        # pkg_check_modules with REQUIRED should already be fatal, but double-check.
        message(FATAL_ERROR "gperftools required for profiling but not found.")
    endif()
endif()

#-------------------------------------------------------
# Spdlog (Logging - Required)
#-------------------------------------------------------
if(MEMILIO_USE_BUNDLED_SPDLOG)
    message(STATUS "Using bundled Spdlog ${MEMILIO_BUNDLED_SPDLOG_VERSION}")
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v${MEMILIO_BUNDLED_SPDLOG_VERSION}
        GIT_SHALLOW    TRUE
        GIT_PROGRESS   TRUE
    )
    # Set FetchContent build options before MakeAvailable
    set(SPDLOG_INSTALL OFF CACHE BOOL "" FORCE)
    set(SPDLOG_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE)
    set(SPDLOG_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(SPDLOG_BUILD_BENCH OFF CACHE BOOL "" FORCE)
    set(SPDLOG_FMT_EXTERNAL OFF CACHE BOOL "" FORCE) # Use spdlog's bundled fmt library
    FetchContent_MakeAvailable(spdlog)

    # Verify the expected target was created by the bundled build
    if(NOT TARGET spdlog::spdlog)
        message(FATAL_ERROR "Bundled spdlog build did not create target spdlog::spdlog.")
    endif()
else()
    message(STATUS "Searching for system Spdlog...")
    find_package(spdlog REQUIRED CONFIG) # Use CMake config mode for modern targets
    message(STATUS "Found system spdlog: ${spdlog_VERSION}")

    # Ensure the modern spdlog::spdlog target exists, creating an alias if needed
    if(TARGET spdlog AND NOT TARGET spdlog::spdlog)
        message(STATUS "Creating alias spdlog::spdlog for system target spdlog.")
        add_library(spdlog::spdlog ALIAS spdlog)
    elseif(NOT TARGET spdlog::spdlog)
         message(FATAL_ERROR "System spdlog found, but required target spdlog::spdlog is missing.")
    endif()
endif()

#-------------------------------------------------------
# Eigen (Linear Algebra - Required)
#-------------------------------------------------------
if(MEMILIO_USE_BUNDLED_EIGEN)
    message(STATUS "Using bundled Eigen ${MEMILIO_BUNDLED_EIGEN_VERSION}")
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        ${MEMILIO_BUNDLED_EIGEN_VERSION}
        GIT_SHALLOW    TRUE
        GIT_PROGRESS   TRUE
    )
    FetchContent_MakeAvailable(eigen) # Populates eigen_SOURCE_DIR

    # Eigen is header-only; manually create the INTERFACE target when fetched
    if(NOT TARGET Eigen3::Eigen)
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
        message(STATUS "Created target Eigen3::Eigen for bundled version.")
    endif()
else()
    # Require system Eigen version >= bundled version for consistency
    set(MEMILIO_MINIMAL_EIGEN_VERSION ${MEMILIO_BUNDLED_EIGEN_VERSION})
    message(STATUS "Searching for system Eigen >= ${MEMILIO_MINIMAL_EIGEN_VERSION}...")
    find_package(Eigen3 ${MEMILIO_MINIMAL_EIGEN_VERSION} REQUIRED NO_MODULE) # Use Eigen3Config.cmake
    message(STATUS "Found system Eigen3: ${EIGEN3_VERSION_STRING}")

    # Verify find_package created the target
    if(NOT TARGET Eigen3::Eigen)
        message(FATAL_ERROR "System Eigen3 found, but required target Eigen3::Eigen is missing.")
    endif()
endif()

#-------------------------------------------------------
# Ipopt (Numerical Optimization - Optional)
#-------------------------------------------------------
if(MEMILIO_ENABLE_IPOPT)
    # Currently, only the bundled version via a specific CMake-enabled fork is supported.
    message(STATUS "Using bundled Ipopt ${MEMILIO_BUNDLED_IPOPT_VERSION} (MEMILIO_ENABLE_IPOPT=ON)")
    FetchContent_Declare(
        Ipopt # Uses a specific CMake-enabled fork
        GIT_REPOSITORY https://git.rwth-aachen.de/avt-svt/public/thirdparty/IpoptCmake.git
        GIT_TAG        ${MEMILIO_BUNDLED_IPOPT_VERSION}
        GIT_SHALLOW    TRUE
        GIT_PROGRESS   TRUE
    )
    FetchContent_MakeAvailable(Ipopt)

    # Verify that the Ipopt target was created successfully
    # Adjust target name checks if the IpoptCmake project uses different ones
    if(TARGET ipopt AND NOT TARGET Ipopt::Ipopt)
        message(STATUS "Bundled Ipopt configured. Creating alias Ipopt::Ipopt.")
        add_library(Ipopt::Ipopt ALIAS ipopt)
    elseif(TARGET Ipopt::Ipopt)
         message(STATUS "Bundled Ipopt configured (target Ipopt::Ipopt found).")
    else()
         # Warn but don't fail by default, allows project to proceed without Ipopt if desired
         message(WARNING "Ipopt target (expected 'ipopt' or 'Ipopt::Ipopt') not found after FetchContent.")
    endif()
endif()

#-------------------------------------------------------
# Boost 
#-------------------------------------------------------
if(MEMILIO_USE_BUNDLED_BOOST)
    message(STATUS "Using bundled Boost ${MEMILIO_BUNDLED_BOOST_VERSION} (building filesystem component)")
    string(REPLACE "." "_" MEMILIO_BOOST_VERSION_UNDERSC "${MEMILIO_BUNDLED_BOOST_VERSION}")
    FetchContent_Declare(
        boost_download # Just downloads the source archive
        URL               https://archives.boost.io/release/${MEMILIO_BUNDLED_BOOST_VERSION}/source/boost_${MEMILIO_BOOST_VERSION_UNDERSC}.tar.gz
        URL_HASH          SHA256=bcce3562880eea4eb6fbf1914a8f3563c21c5956fa066b984dbb3f96c1ca6766
        DOWNLOAD_PROGRESS TRUE
    )
    FetchContent_MakeAvailable(boost_download) # Populates boost_download_SOURCE_DIR

    # Manual build for bundled Boost::filesystem (avoids full Boost build system)
    # 1. Interface target for headers
    add_library(Boost::boost INTERFACE)
    target_include_directories(Boost::boost SYSTEM INTERFACE ${boost_download_SOURCE_DIR})
    # Suppress common Boost header warnings
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(Boost::boost INTERFACE "-Wno-c++20-attribute-extensions")
    endif()

    # 2. Interface target to disable MSVC auto-linking
    add_library(Boost::disable_autolinking INTERFACE)
    if(MSVC)
        target_compile_definitions(Boost::disable_autolinking INTERFACE BOOST_ALL_NO_LIB)
    endif()

    # 3. Static library target for filesystem source files
    add_library(Boost::filesystem STATIC
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/codecvt_error_category.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/directory.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/exception.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/operations.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/path.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/path_traits.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/portability.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/unique_path.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/utf8_codecvt_facet.cpp
        # Optionally add platform-specific sources using generator expressions
        # $<IF:$<PLATFORM_ID:Windows>,${boost_download_SOURCE_DIR}/libs/filesystem/src/windows_file_codecvt.cpp,"">
    )
    target_compile_definitions(Boost::filesystem PUBLIC BOOST_FILESYSTEM_NO_CXX20_ATOMIC_REF)
    if(NOT MSVC)
       target_compile_definitions(Boost::filesystem PUBLIC BOOST_NO_CXX98_FUNCTION_BASE)
    endif()
    target_link_libraries(Boost::filesystem PUBLIC Boost::boost Boost::disable_autolinking)
    set_property(TARGET Boost::filesystem PROPERTY POSITION_INDEPENDENT_CODE ON) # Needed for static libs used in shared libs or executables

    # Signal that Boost (required parts) is available
    set(Boost_FOUND ON CACHE BOOL "Boost found (bundled)" FORCE)
    set(Boost_VERSION_STRING ${MEMILIO_BUNDLED_BOOST_VERSION} CACHE STRING "Boost version (bundled)" FORCE)
    message(STATUS "Bundled Boost::filesystem target created.")

else()
    message(STATUS "Searching for system Boost >= ${MEMILIO_MINIMAL_BOOST_VERSION} (components: filesystem)...")
    find_package(Boost ${MEMILIO_MINIMAL_BOOST_VERSION} REQUIRED COMPONENTS filesystem) # find_package should define Boost::filesystem
    message(STATUS "Found system Boost: ${Boost_VERSION_STRING}")
    # Verify find_package created the target
    if(NOT TARGET Boost::filesystem)
        message(FATAL_ERROR "System Boost found, but required target Boost::filesystem is missing.")
    endif()
endif()

#-------------------------------------------------------
# HDF5 (IO Format - Optional)
#-------------------------------------------------------
# Check MEMILIO_ENABLE_HDF5 option set in the main CMakeLists.txt
set(MEMILIO_HAS_HDF5 FALSE CACHE BOOL "HDF5 support status" FORCE) # Initialize status variable
if(MEMILIO_ENABLE_HDF5)
    message(STATUS "Searching for optional system HDF5 (C library)...")
    # Find quietly using CONFIG mode first, module mode as fallback (standard find_package behavior)
    find_package(HDF5 COMPONENTS C QUIET CONFIG)

    if(HDF5_FOUND AND TARGET HDF5::HDF5) # Check for FOUND status and modern target
        message(STATUS "Found system HDF5: ${HDF5_VERSION} (using target HDF5::HDF5)")
        set(MEMILIO_HAS_HDF5 TRUE CACHE BOOL "HDF5 support status" FORCE)
    else()
        message(STATUS "Optional: System HDF5 not found or target HDF5::HDF5 missing.")
        # MEMILIO_HAS_HDF5 remains FALSE
    endif()
else()
     message(STATUS "HDF5 support explicitly disabled (MEMILIO_ENABLE_HDF5=OFF).")
endif()

#-------------------------------------------------------
# JsonCpp (JSON Parsing - Optional)
#-------------------------------------------------------
# Option MEMILIO_USE_BUNDLED_JSONCPP controls whether to fetch or find
set(MEMILIO_HAS_JSONCPP FALSE CACHE BOOL "JsonCpp support status" FORCE) # Initialize status variable
if(MEMILIO_USE_BUNDLED_JSONCPP)
    message(STATUS "Using bundled JsonCpp ${MEMILIO_BUNDLED_JSONCPP_VERSION}")
    FetchContent_Declare(
        jsoncpp
        URL               https://github.com/open-source-parsers/jsoncpp/archive/${MEMILIO_BUNDLED_JSONCPP_VERSION}.tar.gz
        URL_HASH          SHA256=803bba108191087f8d499b2acfdc7e6e598d45cb93daceaf913248f6ee55335a
        DOWNLOAD_PROGRESS TRUE
    )
    # Set FetchContent build options before MakeAvailable
    set(JSONCPP_WITH_TESTS OFF CACHE BOOL "" FORCE)
    set(JSONCPP_WITH_POST_BUILD_UNITTEST OFF CACHE BOOL "" FORCE)
    set(JSONCPP_INSTALL OFF CACHE BOOL "" FORCE)
    # Match shared/static build type with the main project
    set(BUILD_SHARED_LIBS ${MEMILIO_BUILD_SHARED_LIBS} CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(jsoncpp)

    # Create the standard JsonCpp::JsonCpp alias target if needed
    if(TARGET jsoncpp_static AND NOT TARGET JsonCpp::JsonCpp)
        add_library(JsonCpp::JsonCpp ALIAS jsoncpp_static)
        message(STATUS "Bundled JsonCpp (static) configured. Target JsonCpp::JsonCpp created.")
        set(MEMILIO_HAS_JSONCPP TRUE CACHE BOOL "JsonCpp support status" FORCE)
    elseif(TARGET jsoncpp_lib AND NOT TARGET JsonCpp::JsonCpp)
        add_library(JsonCpp::JsonCpp ALIAS jsoncpp_lib)
        message(STATUS "Bundled JsonCpp (shared) configured. Target JsonCpp::JsonCpp created.")
        set(MEMILIO_HAS_JSONCPP TRUE CACHE BOOL "JsonCpp support status" FORCE)
    elseif(TARGET JsonCpp::JsonCpp) # If subproject already created the namespaced target
         message(STATUS "Bundled JsonCpp configured (Target JsonCpp::JsonCpp exists).")
         set(MEMILIO_HAS_JSONCPP TRUE CACHE BOOL "JsonCpp support status" FORCE)
    else()
        message(WARNING "Failed to find expected target after configuring bundled JsonCpp.")
        # MEMILIO_HAS_JSONCPP remains FALSE
    endif()

else()
    message(STATUS "Searching for optional system JsonCpp...")
    # Find quietly using CONFIG mode first
    find_package(jsoncpp CONFIG QUIET)

    if(TARGET JsonCpp::JsonCpp) # Check for the standard modern target
         set(jsoncpp_version_str "unknown")
         # Try to get version info
         if(jsoncpp_VERSION) set(jsoncpp_version_str ${jsoncpp_VERSION})
         elseif(JsonCpp_VERSION) set(jsoncpp_version_str ${JsonCpp_VERSION}) endif()
         message(STATUS "Found system JsonCpp: ${jsoncpp_version_str} (using target JsonCpp::JsonCpp)")
         set(MEMILIO_HAS_JSONCPP TRUE CACHE BOOL "JsonCpp support status" FORCE)
    else()
         message(STATUS "Optional: System JsonCpp not found or JsonCpp::JsonCpp target missing.")
         # MEMILIO_HAS_JSONCPP remains FALSE
    endif()
endif()

#-------------------------------------------------------
# Random123 (RNG - Required, Header-Only)
#-------------------------------------------------------
# Always use the bundled version for this simple header-only library
message(STATUS "Using bundled Random123 ${MEMILIO_BUNDLED_RANDOM123_VERSION}")
FetchContent_Declare(
    Random123
    GIT_REPOSITORY https://github.com/DEShawResearch/random123
    GIT_TAG        ${MEMILIO_BUNDLED_RANDOM123_VERSION}
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
)
FetchContent_MakeAvailable(Random123) # Populates random123_SOURCE_DIR

# Create the INTERFACE target
if(NOT TARGET Random123::Random123)
    add_library(Random123::Random123 INTERFACE)
    target_include_directories(Random123::Random123 SYSTEM INTERFACE ${random123_SOURCE_DIR}/include)
    message(STATUS "Bundled Random123 configured. Target Random123::Random123 created.")
endif()
