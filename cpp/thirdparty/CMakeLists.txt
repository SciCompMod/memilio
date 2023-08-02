# Versions of the bundled libraries
# If you like to upgrade, just change the number
set(MEMILIO_EIGEN_VERSION "3.3.9")
set(MEMILIO_SPDLOG_VERSION "1.11.0")
set(MEMILIO_JSONCPP_VERSION "1.9.5")

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
    target_include_directories(eigen INTERFACE ${eigen_SOURCE_DIR})
    add_library(Eigen3::Eigen ALIAS eigen)
else()
    find_package(Eigen3 ${MEMILIO_EIGEN_VERSION} REQUIRED NO_MODULE)
endif()

# ## BOOST
if(MEMILIO_USE_BUNDLED_BOOST)
    include(BuildBoost)
else()
    find_package(Boost REQUIRED COMPONENTS outcome optional filesystem)
endif(MEMILIO_USE_BUNDLED_BOOST)

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

if (MEMILIO_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()
