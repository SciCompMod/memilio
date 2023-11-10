# Versions of the bundled libraries
# If you like to upgrade, just change the number
set(MEMILIO_SPDLOG_VERSION "1.11.0")
set(MEMILIO_EIGEN_VERSION "3.3.9")
set(MEMILIO_JSONCPP_VERSION "1.9.5")
set(MEMILIO_BOOST_VERSION "1.75.0")
include(FetchContent)
#TODO: There can be a heavy reduction in setup if we set some things right here
set(MEMILIO_RANDOM123_VERSION "v1.14.0")

# ## SPDLOG
MESSAGE(STATUS "Searching for SPDLOG and downloading if not found")
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v${MEMILIO_SPDLOG_VERSION}
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    FIND_PACKAGE_ARGS ${MEMILIO_SPDLOG_VERSION}
)

set(SPDLOG_INSTALL ON)

FetchContent_MakeAvailable(spdlog)


# ## EIGEN
MESSAGE(STATUS "Searching for Eigen and downloading if not found")
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG ${MEMILIO_EIGEN_VERSION}
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
    # FIND_PACKAGE_ARGS ${MEMILIO_EIGEN_VERSION} NAMES eigen NO_MODULE
)

FetchContent_MakeAvailable(eigen)
add_library(Eigen3::Eigen ALIAS eigen)

# ## BOOST
MESSAGE(STATUS "Searching for Boost and downloading if not found")
include(CMake/thirdparty/BuildBoost)

# ## HDF5
MESSAGE(STATUS "Looking for HDF5")
find_package(HDF5 COMPONENTS C)

if(HDF5_FOUND)
    set(MEMILIO_HAS_HDF5 ON)
else()
    message(WARNING "HDF5 was not found. Memilio will be built without some IO features. Install HDF5 Libraries and set the HDF5_DIR cmake variable to the directory containing the hdf5-config.cmake file to build with HDF5.")
endif()

# ## JSONCPP
MESSAGE(STATUS "Searching for Jsoncpp and downloading if not found")


include(FetchContent)
FetchContent_Declare(
    jsoncpp
    URL https://github.com/open-source-parsers/jsoncpp/archive/${MEMILIO_JSONCPP_VERSION}.tar.gz
)

set(JSONCPP_WITH_TESTS "Compile and (for jsoncpp_check) run JsonCpp test executables" OFF)
set(JSONCPP_WITH_POST_BUILD_UNITTEST "Automatically run unit-tests as a post build step" ON)
FetchContent_MakeAvailable(jsoncpp)
# unset global cache variables to avoid clashes with our code
unset(BUILD_OBJECT_LIBS CACHE)
unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY CACHE)
unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY CACHE)
unset(CMAKE_PDB_OUTPUT_DIRECTORY CACHE)
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY CACHE)
if(BUILD_SHARED_LIBS)
    add_library(JsonCpp::JsonCpp ALIAS jsoncpp_lib)
else()
    add_library(JsonCpp::JsonCpp ALIAS jsoncpp_static)
endif()
if(TARGET JsonCpp::JsonCpp)
    set(MEMILIO_HAS_JSONCPP ON)
else()
    message(WARNING "JsonCpp was not found. Memilio will be built without some IO features. 
        Set CMake variable MEMILIO_USE_BUNDLED_JSONCPP to ON or install JsonCpp and set the jsoncpp_DIR cmake variable 
        to the directory containing the jsoncppConfig.cmake file to build with JsonCpp.")
endif()


## MPI
if (MEMILIO_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
endif()

if (MEMILIO_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()

#Random123 library for random number generators
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
