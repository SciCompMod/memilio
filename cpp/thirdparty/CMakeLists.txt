# Versions of the bundled libraries
# If you like to upgrade, just change the number
set(MEMILIO_EIGEN_VERSION "3.4.0")
set(MEMILIO_SPDLOG_VERSION "1.17.0")
set(MEMILIO_BOOST_VERSION "1.84.0")
set(MEMILIO_MINIMAL_BOOST_VERSION "1.76.0")
set(MEMILIO_JSONCPP_VERSION "1.9.6")
set(MEMILIO_RANDOM123_VERSION "v1.14.0")
set(MEMILIO_IPOPT_VERSION "3.14.12")
set(MEMILIO_SBML_VERSION "5.20.2")

include(FetchContent)

# Bundled header-only libraries are wrapped in our own INTERFACE targets. When installing, their headers are installed
# to this directory and the targets are exported together with memilio.
set(MEMILIO_THIRDPARTY_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/memilio/thirdparty")
set(MEMILIO_BUNDLED_TARGETS)

# Dependencies are EXCLUDE_FROM_ALL, so only the targets that memilio uses are built. Install rules of such dependencies
# are skipped, so spdlog and jsoncpp, which are installed next to memilio, only use it if memilio is not installed.
if(MEMILIO_ENABLE_INSTALL)
    set(MEMILIO_EXCLUDE_FROM_ALL_IF_NOT_INSTALLED)
else()
    set(MEMILIO_EXCLUDE_FROM_ALL_IF_NOT_INSTALLED EXCLUDE_FROM_ALL)
endif()

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
if(MEMILIO_USE_BUNDLED_SPDLOG)
    message(STATUS "Downloading Spdlog library")

    set(SPDLOG_INSTALL ${MEMILIO_ENABLE_INSTALL})
    FetchContent_Declare(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v${MEMILIO_SPDLOG_VERSION}
        GIT_SHALLOW TRUE
        SYSTEM
        ${MEMILIO_EXCLUDE_FROM_ALL_IF_NOT_INSTALLED}
    )
    FetchContent_MakeAvailable(spdlog)
else()
    find_package(spdlog REQUIRED)
endif()

# ## EIGEN
if(MEMILIO_USE_BUNDLED_EIGEN)
    message(STATUS "Downloading Eigen library")

    # Set SOURCE_SUBDIR to prevent MakeAvailable from calling add_subdirectory. We do not want to run Eigen's
    # CMakeLists.txt, as it is not set up to be used like that. Instead, we need to manually add a library with the
    # source files below (this is sufficient, as Eigen is a header-only library)
    FetchContent_Declare(eigen
        SOURCE_SUBDIR this_directory_name_does_not_exist
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG ${MEMILIO_EIGEN_VERSION}
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(eigen)

    add_library(eigen INTERFACE)
    target_include_directories(eigen SYSTEM INTERFACE
        $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${MEMILIO_THIRDPARTY_INCLUDEDIR}>
    )
    add_library(Eigen3::Eigen ALIAS eigen)
    list(APPEND MEMILIO_BUNDLED_TARGETS eigen)
    if(MEMILIO_ENABLE_INSTALL)
        install(DIRECTORY ${eigen_SOURCE_DIR}/Eigen ${eigen_SOURCE_DIR}/unsupported
            DESTINATION ${MEMILIO_THIRDPARTY_INCLUDEDIR})
    endif()
else()
    find_package(Eigen3 ${MEMILIO_EIGEN_VERSION} REQUIRED NO_MODULE)
endif()

# ## IPOPT (numerial optimization)
if(MEMILIO_ENABLE_IPOPT)
    message(STATUS "Downloading Ipopt library")

    FetchContent_Declare(Ipopt
        GIT_REPOSITORY https://git.rwth-aachen.de/avt-svt/public/thirdparty/IpoptCmake.git
        GIT_TAG ${MEMILIO_IPOPT_VERSION}
        GIT_SHALLOW TRUE
        SYSTEM
        EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(Ipopt)
endif()

# ## BOOST
if(MEMILIO_USE_BUNDLED_BOOST)
    message(STATUS "Downloading Boost library")

    string(REPLACE "." "_" MEMILIO_BOOST_VERSION_UNDERSC "${MEMILIO_BOOST_VERSION}")

    FetchContent_Declare(boost
        # don't use the URL from github, that download isn't complete and requires more setup (subrepositories, bootstrapping)
        URL https://archives.boost.io/release/${MEMILIO_BOOST_VERSION}/source/boost_${MEMILIO_BOOST_VERSION_UNDERSC}.tar.gz
    )
    FetchContent_MakeAvailable(boost)

    add_library(boost INTERFACE)
    add_library(Boost::boost ALIAS boost)
    target_include_directories(boost SYSTEM INTERFACE
        $<BUILD_INTERFACE:${boost_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${MEMILIO_THIRDPARTY_INCLUDEDIR}>
    )

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

    list(APPEND MEMILIO_BUNDLED_TARGETS boost boost_disable_autolink)
    if(MEMILIO_ENABLE_INSTALL)
        install(DIRECTORY ${boost_SOURCE_DIR}/boost DESTINATION ${MEMILIO_THIRDPARTY_INCLUDEDIR})
    endif()

    set(Boost_LIBRARIES Boost::boost)
    set(Boost_FOUND ON)
else()
    find_package(Boost ${MEMILIO_MINIMAL_BOOST_VERSION}...${MEMILIO_BOOST_VERSION} REQUIRED COMPONENTS outcome optional)
endif()

# ## HDF5
if(MEMILIO_ENABLE_HDF5)
    find_package(HDF5 COMPONENTS C)

    if(HDF5_FOUND)
        set(MEMILIO_HAS_HDF5 ON)
    else()
        message(WARNING "HDF5 was not found. Memilio will be built without some IO features. Install HDF5 Libraries and set the HDF5_DIR cmake variable to the directory containing the hdf5-config.cmake file to build with HDF5.")
    endif()
endif()

# ## JSONCPP
if(MEMILIO_USE_BUNDLED_JSONCPP)
    message(STATUS "Downloading JsonCpp library")

    # set jsoncpp configurations, only build the library variant that memilio links against
    set(JSONCPP_WITH_TESTS OFF)
    set(JSONCPP_WITH_POST_BUILD_UNITTEST OFF)
    set(JSONCPP_WITH_CMAKE_PACKAGE ${MEMILIO_ENABLE_INSTALL})
    set(BUILD_OBJECT_LIBS OFF)
    if(BUILD_SHARED_LIBS)
        set(BUILD_STATIC_LIBS OFF)
    else()
        set(BUILD_STATIC_LIBS ON)
    endif()

    FetchContent_Declare(jsoncpp
        URL https://github.com/open-source-parsers/jsoncpp/archive/${MEMILIO_JSONCPP_VERSION}.tar.gz
        SYSTEM
        ${MEMILIO_EXCLUDE_FROM_ALL_IF_NOT_INSTALLED}
    )
    FetchContent_MakeAvailable(jsoncpp)

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


if(MEMILIO_ENABLE_SBML)
    find_package(sbml REQUIRED ${MEMILIO_SBML_VERSION} CONFIG)
    message(STATUS "Found libSBML: ${libsbml_INCLUDE_DIRS}")
endif()

if(MEMILIO_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    message(STATUS "Enabled MPI ${MPI_CXX_VERSION}")
endif()

if(MEMILIO_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    message(STATUS "Enabled OpenMP ${OpenMP_CXX_VERSION}")
endif()

# Random123 library for random number generators
message(STATUS "Downloading Random123 library")

FetchContent_Declare(Random123
    GIT_REPOSITORY https://github.com/DEShawResearch/random123
    GIT_TAG ${MEMILIO_RANDOM123_VERSION}
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Random123)

add_library(Random123 INTERFACE)
target_include_directories(Random123 SYSTEM INTERFACE
    $<BUILD_INTERFACE:${random123_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${MEMILIO_THIRDPARTY_INCLUDEDIR}>
)
list(APPEND MEMILIO_BUNDLED_TARGETS Random123)
if(MEMILIO_ENABLE_INSTALL)
    install(DIRECTORY ${random123_SOURCE_DIR}/include/Random123 DESTINATION ${MEMILIO_THIRDPARTY_INCLUDEDIR})
endif()
