# Versions of the bundled libraries
# If you like to upgrade, just change the number
set(MEMILIO_EIGEN_VERSION "3.4.0")
set(MEMILIO_SPDLOG_VERSION "1.15.0")
set(MEMILIO_BOOST_VERSION "1.84.0")
set(MEMILIO_MINIMAL_BOOST_VERSION "1.76.0")
set(MEMILIO_JSONCPP_VERSION "1.9.6")
set(MEMILIO_RANDOM123_VERSION "v1.14.0")
set(MEMILIO_IPOPT_VERSION "3.14.12")

# Gperftools for profiling; must be first, so that libraries are included in the profile
if(MEMILIO_ENABLE_PROFILING)
    if(NOT UNIX)
        message(FATAL_ERROR "Profiling with gperftools currently only supported on unix systems.")
    endif()

    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GPERFTOOLS REQUIRED libprofiler)
    link_libraries(${GPERFTOOLS_LINK_LIBRARIES}) # link globally so all libs are included
    include_directories(${GPERFTOOLS_INCLUDE_DIRS})
endif()

# ## SPDLOG
set(SPDLOG_INSTALL ON)

if(MEMILIO_USE_BUNDLED_SPDLOG)
    message(STATUS "Downloading Spdlog library")

    include(FetchContent)
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v${MEMILIO_SPDLOG_VERSION}
    )
    FetchContent_GetProperties(spdlog)

    if(NOT spdlog_POPULATED)
        FetchContent_Populate(spdlog)
        add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()

else()
    find_package(spdlog REQUIRED)
endif()

# ## EIGEN
if(MEMILIO_USE_BUNDLED_EIGEN)
    message(STATUS "Downloading Eigen library")

    include(FetchContent)
    FetchContent_Declare(eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG ${MEMILIO_EIGEN_VERSION})
    FetchContent_GetProperties(eigen)

    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
    endif()

    add_library(eigen INTERFACE)
    target_include_directories(eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
    add_library(Eigen3::Eigen ALIAS eigen)
else()
    find_package(Eigen3 ${MEMILIO_EIGEN_VERSION} REQUIRED NO_MODULE)
endif()

# ## IPOPT (numerial optimization)
if(MEMILIO_ENABLE_IPOPT)
    message(STATUS "Downloading Ipopt library")

    include(FetchContent)
    FetchContent_Declare(Ipopt
        GIT_REPOSITORY https://git.rwth-aachen.de/avt-svt/public/thirdparty/IpoptCmake.git
        GIT_TAG ${MEMILIO_IPOPT_VERSION})

    FetchContent_MakeAvailable(Ipopt)
endif()

# ## BOOST
if(MEMILIO_USE_BUNDLED_BOOST)
    message(STATUS "Downloading Boost library")

    string(REPLACE "." "_" MEMILIO_BOOST_VERSION_UNDERSC "${MEMILIO_BOOST_VERSION}")

    include(FetchContent)
    FetchContent_Declare(boost

        # don't use the URL from github, that download isn't complete and requires more setup (subrepositories, bootstrapping)
        URL https://archives.boost.io/release/${MEMILIO_BOOST_VERSION}/source/boost_${MEMILIO_BOOST_VERSION_UNDERSC}.tar.gz
    )
    FetchContent_GetProperties(boost)

    if(NOT boost_POPULATED)
        FetchContent_Populate(boost)
    endif()

    add_library(boost INTERFACE)
    add_dependencies(boost boost-bootstrap)
    add_library(Boost::boost ALIAS boost)
    target_include_directories(boost SYSTEM INTERFACE $<BUILD_INTERFACE:${boost_SOURCE_DIR}>)

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
            target_compile_options(boost INTERFACE "-Wno-c++20-attribute-extensions")
        else()
            target_compile_options(boost INTERFACE "-Wno-c++20-extensions")
        endif()
    endif()

    add_library(boost_disable_autolink INTERFACE)
    target_compile_definitions(boost_disable_autolink INTERFACE BOOST_ALL_NO_LIB)
    add_library(Boost::disable_autolinking ALIAS boost_disable_autolink)

    add_library(boost_filesystem STATIC
        ${boost_SOURCE_DIR}/libs/filesystem/src/codecvt_error_category.cpp
        ${boost_SOURCE_DIR}/libs/filesystem/src/directory.cpp
        ${boost_SOURCE_DIR}/libs/filesystem/src/exception.cpp
        ${boost_SOURCE_DIR}/libs/filesystem/src/operations.cpp
        ${boost_SOURCE_DIR}/libs/filesystem/src/path.cpp
        ${boost_SOURCE_DIR}/libs/filesystem/src/path_traits.cpp
        ${boost_SOURCE_DIR}/libs/filesystem/src/portability.cpp
        ${boost_SOURCE_DIR}/libs/filesystem/src/unique_path.cpp
        ${boost_SOURCE_DIR}/libs/filesystem/src/utf8_codecvt_facet.cpp
        ${boost_SOURCE_DIR}/libs/filesystem/src/windows_file_codecvt.cpp
    )

    # Ensure that the boost atomic library is used instead of the standard atomic library, where some functionality is only available as of C++20.
    target_compile_definitions(boost_filesystem PUBLIC BOOST_FILESYSTEM_NO_CXX20_ATOMIC_REF)

    target_link_libraries(boost_filesystem PUBLIC boost_disable_autolink boost)
    set_property(TARGET boost_filesystem PROPERTY POSITION_INDEPENDENT_CODE ON)
    add_library(Boost::filesystem ALIAS boost_filesystem)

    if(NOT MSVC) # on gcc and apple clang we need to define BOOST_NO_CXX98_FUNCTION_BASE because a deprecated function is sometimes used in boost
        target_compile_definitions(boost_filesystem PUBLIC BOOST_NO_CXX98_FUNCTION_BASE)
    endif()

    set(Boost_LIBRARIES Boost::boost Boost::filesystem)
    set(Boost_FOUND ON)

else()
    find_package(Boost ${MEMILIO_MINIMAL_BOOST_VERSION}...${MEMILIO_BOOST_VERSION} REQUIRED COMPONENTS outcome optional filesystem)
endif()

# ## HDF5
find_package(HDF5 COMPONENTS C)

if(HDF5_FOUND)
    set(MEMILIO_HAS_HDF5 ON)
else()
    message(WARNING "HDF5 was not found. Memilio will be built without some IO features. Install HDF5 Libraries and set the HDF5_DIR cmake variable to the directory containing the hdf5-config.cmake file to build with HDF5.")
endif()

# ## JSONCPP
if(MEMILIO_USE_BUNDLED_JSONCPP)
    message(STATUS "Downloading jsoncpp library")

    include(FetchContent)
    FetchContent_Declare(
        jsoncpp
        URL https://github.com/open-source-parsers/jsoncpp/archive/${MEMILIO_JSONCPP_VERSION}.tar.gz
    )
    FetchContent_GetProperties(jsoncpp)

    if(NOT jsoncpp_POPULATED)
        FetchContent_Populate(jsoncpp)

        # set jsoncpp configurations
        set(JSONCPP_WITH_TESTS "Compile and (for jsoncpp_check) run JsonCpp test executables" OFF)
        set(JSONCPP_WITH_POST_BUILD_UNITTEST "Automatically run unit-tests as a post build step" ON)

        add_subdirectory(${jsoncpp_SOURCE_DIR} ${jsoncpp_BINARY_DIR} EXCLUDE_FROM_ALL)

        # unset global cache variables to avoid clashes with our code
        unset(BUILD_OBJECT_LIBS CACHE)
        unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY CACHE)
        unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY CACHE)
        unset(CMAKE_PDB_OUTPUT_DIRECTORY CACHE)
        unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY CACHE)
    endif()

    if(BUILD_SHARED_LIBS)
        add_library(JsonCpp::JsonCpp ALIAS jsoncpp_lib)
    else()
        add_library(JsonCpp::JsonCpp ALIAS jsoncpp_static)
    endif()
else()
    find_package(jsoncpp CONFIG)
endif()

if(TARGET JsonCpp::JsonCpp)
    set(MEMILIO_HAS_JSONCPP ON)
else()
    message(WARNING "JsonCpp was not found. Memilio will be built without some IO features. 
        Set CMake variable MEMILIO_USE_BUNDLED_JSONCPP to ON or install JsonCpp and set the jsoncpp_DIR cmake variable 
        to the directory containing the jsoncppConfig.cmake file to build with JsonCpp.")
endif()

if(MEMILIO_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()

if(MEMILIO_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

# Random123 library for random number generators
message(STATUS "Downloading Random123 library")

include(FetchContent)
FetchContent_Declare(Random123
    GIT_REPOSITORY https://github.com/DEShawResearch/random123
    GIT_TAG ${MEMILIO_RANDOM123_VERSION}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND "")
FetchContent_GetProperties(Random123)

if(NOT Random123_POPULATED)
    FetchContent_Populate(Random123)
endif()

add_library(Random123 INTERFACE)
target_include_directories(Random123 INTERFACE ${random123_SOURCE_DIR}/include)
