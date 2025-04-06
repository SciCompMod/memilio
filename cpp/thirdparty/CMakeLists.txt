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

    if(GPERFTOOLS_FOUND AND NOT TARGET memilio::gperftools)
        add_library(memilio::gperftools INTERFACE IMPORTED)
        target_include_directories(memilio::gperftools SYSTEM INTERFACE ${GPERFTOOLS_INCLUDE_DIRS}) 
        target_link_libraries(memilio::gperftools INTERFACE ${GPERFTOOLS_LINK_LIBRARIES})
        message(STATUS "Profiling enabled: Found system gperftools. Created target memilio::gperftools.")
    elseif(NOT GPERFTOOLS_FOUND)
         message(FATAL_ERROR "gperftools required for profiling but not found.")
    endif()
endif()

#-------------------------------------------------------
# Spdlog (Logging - Required)
#-------------------------------------------------------
set(MEMILIO_HAS_SPDLOG FALSE CACHE BOOL "Spdlog support status" FORCE) # Initialize status variable
if(MEMILIO_USE_BUNDLED_SPDLOG)
    message(STATUS "Using bundled Spdlog ${MEMILIO_BUNDLED_SPDLOG_VERSION}")
    FetchContent_Declare(
        spdlog
        QUIET
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
        if(TARGET spdlog)
            message(STATUS "Creating alias spdlog::spdlog for bundled target spdlog.")
            add_library(spdlog::spdlog ALIAS spdlog)
            set(MEMILIO_HAS_SPDLOG TRUE CACHE BOOL "Spdlog support status" FORCE)
        else()
            message(FATAL_ERROR "Bundled spdlog build did not create expected target spdlog or spdlog::spdlog.")
        endif()
    else()
        set(MEMILIO_HAS_SPDLOG TRUE CACHE BOOL "Spdlog support status" FORCE) # Update status variable
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
    set(MEMILIO_HAS_SPDLOG TRUE CACHE BOOL "Spdlog support status" FORCE) # Update status variable
endif()

#-------------------------------------------------------
# Eigen (Linear Algebra - Required)
#-------------------------------------------------------
set(MEMILIO_HAS_EIGEN FALSE CACHE BOOL "EIGEN support status" FORCE)
if(MEMILIO_USE_BUNDLED_EIGEN)
    message(STATUS "Using bundled Eigen ${MEMILIO_BUNDLED_EIGEN_VERSION}")
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        ${MEMILIO_BUNDLED_EIGEN_VERSION}
        GIT_SHALLOW    TRUE
        GIT_PROGRESS   TRUE
        QUIET
    )
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE) 
    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE) 
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)     
    set(EIGEN_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen) # Populates eigen_SOURCE_DIR
    if(NOT TARGET Eigen3::Eigen)
        add_library(Eigen3::Eigen INTERFACE)
        target_include_directories(Eigen3::Eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
        message(STATUS "Created target Eigen3::Eigen for bundled version.")
    endif()
    set(MEMILIO_HAS_EIGEN TRUE CACHE BOOL "EIGEN support status" FORCE)
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
    set(MEMILIO_HAS_EIGEN TRUE CACHE BOOL "EIGEN support status" FORCE)
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
        set(IPOPT_FOUND TRUE CACHE BOOL "Ipopt found (bundled)") 
    elseif(TARGET Ipopt::Ipopt)
         message(STATUS "Bundled Ipopt configured (target Ipopt::Ipopt found).")
         set(IPOPT_FOUND TRUE CACHE BOOL "Ipopt found (bundled)") 
    else()
         # Warn but don't fail by default, allows project to proceed without Ipopt if desired
         message(WARNING "Ipopt target (expected 'ipopt' or 'Ipopt::Ipopt') not found after FetchContent.")
    endif()
endif()

#-------------------------------------------------------
# Boost 
#-------------------------------------------------------
set(BOOST_FOUND FALSE CACHE BOOL "Boost found status" FORCE)
if(MEMILIO_USE_BUNDLED_BOOST)
    message(STATUS "Using bundled Boost ${MEMILIO_BUNDLED_BOOST_VERSION} (building filesystem component)")
    string(REPLACE "." "_" MEMILIO_BOOST_VERSION_UNDERSC "${MEMILIO_BUNDLED_BOOST_VERSION}")
    FetchContent_Declare(
        boost_download 
        DOWNLOAD_PROGRESS TRUE
        URL               https://archives.boost.io/release/${MEMILIO_BUNDLED_BOOST_VERSION}/source/boost_${MEMILIO_BOOST_VERSION_UNDERSC}.tar.gz
        URL_HASH          SHA256=a5800f405508f5df8114558ca9855d2640a2de8f0445f051fa1c7c3383045724
        
    )
    FetchContent_MakeAvailable(boost_download) 

    #NOTE: Manually building Boost::filesystem to avoid full Boost build system dependency.
    # 1. Interface target for headers
    add_library(boost INTERFACE)
    target_include_directories(boost SYSTEM INTERFACE ${boost_download_SOURCE_DIR})
    add_library(Boost::boost ALIAS boost)
    # Suppress common Boost header warnings
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(boost INTERFACE "-Wno-c++20-attribute-extensions")
    endif()

    # 2. Interface target to disable MSVC auto-linking
    add_library(boost_disable_autolink INTERFACE)
    add_library(Boost::disable_autolinking ALIAS boost_disable_autolink)
    if(MSVC)
        target_compile_definitions(boost_disable_autolink INTERFACE BOOST_ALL_NO_LIB)
    endif()

    # 3. Static library target for filesystem source files
    add_library(boost_filesystem STATIC
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/codecvt_error_category.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/directory.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/exception.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/operations.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/path.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/path_traits.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/portability.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/unique_path.cpp
        ${boost_download_SOURCE_DIR}/libs/filesystem/src/utf8_codecvt_facet.cpp
        $<IF:$<PLATFORM_ID:Windows>,${boost_download_SOURCE_DIR}/libs/filesystem/src/windows_file_codecvt.cpp,>
    )
    target_compile_definitions(boost_filesystem PUBLIC BOOST_FILESYSTEM_NO_CXX20_ATOMIC_REF)
    if(NOT MSVC)
       target_compile_definitions(boost_filesystem PUBLIC BOOST_NO_CXX98_FUNCTION_BASE)
    endif()
    target_link_libraries(boost_filesystem PUBLIC boost boost_disable_autolink)
    set_property(TARGET boost_filesystem PROPERTY POSITION_INDEPENDENT_CODE ON) # Needed for static libs used in shared libs or executables
    add_library(Boost::filesystem ALIAS boost_filesystem)

    # Signal that Boost (required parts) is available
    set(BOOST_FOUND TRUE CACHE BOOL "Boost found (bundled)" FORCE)
    set(Boost_VERSION_STRING ${MEMILIO_BUNDLED_BOOST_VERSION} CACHE STRING "Boost version (bundled)" FORCE)
    message(VERBOSE "Bundled Boost::filesystem target created.")

else()
    message(STATUS "Searching for system Boost >= ${MEMILIO_MINIMAL_BOOST_VERSION} (components: filesystem)...")
    find_package(Boost ${MEMILIO_MINIMAL_BOOST_VERSION} REQUIRED COMPONENTS filesystem) # find_package should define Boost::filesystem
    message(STATUS "Found system Boost: ${Boost_VERSION_STRING}")
    # Verify find_package created the target
    if(TARGET Boost::filesystem)
        message(FATAL_ERROR "System Boost found, but required target Boost::filesystem is missing.")
        set(BOOST_FOUND TRUE CACHE BOOL "Boost found (system)" FORCE) 
    endif()
endif()

#-------------------------------------------------------
# HDF5 (IO Format)
#-------------------------------------------------------
# Check MEMILIO_ENABLE_HDF5 option set in the main CMakeLists.txt
set(MEMILIO_HAS_HDF5 FALSE CACHE BOOL "HDF5 support status" FORCE) # Initialize status variable
message(STATUS "Searching for optional system HDF5 (C library)...")
find_package(HDF5 COMPONENTS C QUIET)
if(HDF5_FOUND)
    if (TARGET HDF5::HDF5)
        message(STATUS "Found system HDF5: ${HDF5_VERSION} (using target HDF5::HDF5)")
        set(MEMILIO_HAS_HDF5 TRUE CACHE BOOL "HDF5 support status" FORCE)
    elseif(TARGET hdf5) # Check for common non-namespaced target from older find modules
        message(STATUS "Found system HDF5: ${HDF5_VERSION} (target hdf5). Creating alias HDF5::HDF5.")
        add_library(HDF5::HDF5 ALIAS hdf5)
        set(MEMILIO_HAS_HDF5 TRUE CACHE BOOL "HDF5 support status" FORCE)
    else()
         message(WARNING "System HDF5 found, but required target (HDF5::HDF5 or hdf5) missing.")
    endif()
else()
    message(WARNING "System HDF5 not found.")
endif()

#-------------------------------------------------------
# JsonCpp (JSON Parsing)
#-------------------------------------------------------
# Option MEMILIO_USE_BUNDLED_JSONCPP controls whether to fetch or find
set(MEMILIO_HAS_JSONCPP FALSE CACHE BOOL "JsonCpp support status" FORCE) # Initialize status variable
if(MEMILIO_USE_BUNDLED_JSONCPP)
    message(STATUS "Using bundled JsonCpp ${MEMILIO_BUNDLED_JSONCPP_VERSION}")
    FetchContent_Declare(
        jsoncpp
        QUIET
        DOWNLOAD_PROGRESS TRUE
        URL               https://github.com/open-source-parsers/jsoncpp/archive/${MEMILIO_BUNDLED_JSONCPP_VERSION}.tar.gz
        URL_HASH          SHA256=f93b6dd7ce796b13d02c108bc9f79812245a82e577581c4c9aabe57075c90ea2
        
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
    endif()

else()
    message(STATUS "Searching for optional system JsonCpp...")
    # Find quietly using CONFIG mode first
    find_package(jsoncpp CONFIG QUIET)

    if(TARGET JsonCpp::JsonCpp) # Check for the standard modern target
         set(jsoncpp_version_str "unknown")
         if(jsoncpp_VERSION) 
            set(jsoncpp_version_str ${jsoncpp_VERSION})
         elseif(JsonCpp_VERSION) 
            set(jsoncpp_version_str ${JsonCpp_VERSION}) 
        endif()
         message(STATUS "Found system JsonCpp: ${jsoncpp_version_str} (using target JsonCpp::JsonCpp)")
         set(MEMILIO_HAS_JSONCPP TRUE CACHE BOOL "JsonCpp support status" FORCE)
    else()
         message(STATUS "Optional: System JsonCpp not found or JsonCpp::JsonCpp target missing.")
    endif()
endif()

#-------------------------------------------------------
# Random123
#-------------------------------------------------------
# Always use the bundled version for this simple header-only library
message(STATUS "Using bundled Random123 ${MEMILIO_BUNDLED_RANDOM123_VERSION}")
FetchContent_Declare(
    Random123
    QUIET
    GIT_REPOSITORY https://github.com/DEShawResearch/random123
    GIT_TAG        ${MEMILIO_BUNDLED_RANDOM123_VERSION}
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
    SOURCE_SUBDIR ""
)
FetchContent_MakeAvailable(Random123) 

# Create the INTERFACE target
if(NOT TARGET Random123)
    add_library(Random123 INTERFACE)
    message(STATUS ${FETCHCONTENT_BASE_DIR})
    target_include_directories(Random123 SYSTEM INTERFACE ${random123_SOURCE_DIR}/include)
    message(VERBOSE "Bundled Random123 configured. Target Random123 created.")  
endif()


