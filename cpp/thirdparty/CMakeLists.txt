option(MEMILIO_USE_CONAN "Use Conan package manager to resolve dependencies." ON)
option(MEMILIO_USE_SYSTEM_BOOST "Use Boost installed on the system instead of Conan package." OFF)
option(MEMILIO_USE_SYSTEM_EIGEN "Use Eigen3 installed on the system instead of Conan package." OFF)
option(MEMILIO_USE_SYSTEM_SPDLOG "Use spdlog installed on the system instead of Conan package." OFF)
option(MEMILIO_USE_SYSTEM_JSONCPP "Use JsonCpp installed on the system instead of Conan package." OFF)
option(MEMILIO_USE_SYSTEM_HDF5 "Use HDF5 installed on the system instead of Conan package." OFF)

set(_MEMILIO_BOOST_VERSION 1.75.0)
set(_MEMILIO_EIGEN_VERSION 3.3.9)
set(_MEMILIO_SPDLOG_VERSION 1.5.0)
set(_MEMILIO_JSONCPP_VERSION 1.9.5)
set(_MEMILIO_HDF5_VERSION 1.12.0)
set(_MEMILIO_GTEST_VERSION 1.10.0)

if (MEMILIO_USE_CONAN)
    message(STATUS "Memilio: Checking required dependencies.")
    include(conan)

    set(_MEMILIO_REQUIRES "")
    if (NOT MEMILIO_USE_SYSTEM_BOOST)
        list(APPEND _MEMILIO_REQUIRES boost/${_MEMILIO_BOOST_VERSION})
    endif()
    if (NOT MEMILIO_USE_SYSTEM_EIGEN)
        list(APPEND _MEMILIO_REQUIRES eigen/${_MEMILIO_EIGEN_VERSION})
    endif()
    if (NOT MEMILIO_USE_SYSTEM_SPDLOG)
        list(APPEND _MEMILIO_REQUIRES spdlog/${_MEMILIO_SPDLOG_VERSION})
    endif()
    if (NOT MEMILIO_USE_SYSTEM_JSONCPP)
        list(APPEND _MEMILIO_REQUIRES jsoncpp/${_MEMILIO_JSONCPP_VERSION})
    endif()
    if (NOT MEMILIO_USE_SYSTEM_HDF5)
        list(APPEND _MEMILIO_REQUIRES hdf5/${_MEMILIO_HDF5_VERSION})
    endif()
    if (MEMILIO_BUILD_TESTS)
        list(APPEND _MEMILIO_REQUIRES gtest/${_MEMILIO_GTEST_VERSION})
    endif()

    # configure conan
    conan_cmake_configure(
        REQUIRES ${_MEMILIO_REQUIRES}
        GENERATORS cmake_find_package_multi)

    # install conan packages
    get_property(_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if (_MULTI_CONFIG)
        message(STATUS "Memilio: Downloading packages for all configurations.")
        foreach(TYPE ${CMAKE_CONFIGURATION_TYPES})
            conan_cmake_autodetect(_MEMILIO_SETTINGS BUILD_TYPE ${TYPE})
            conan_cmake_install(
                PATH_OR_REFERENCE .
                BUILD missing
                REMOTE conancenter
                SETTINGS ${_MEMILIO_SETTINGS})
        endforeach()
    else()
        message(STATUS "Memilio: Downloading packages.")
        if (NOT CMAKE_BUILD_TYPE)
            message(STATUS "MEMILIO WORKAROUND: Conan requires CMAKE_BUILD_TYPE to be set. Setting default to Release. Overwrite with e.g. -DCMAKE_BUILD_TYPE=<Debug>")
            set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
        endif()
        conan_cmake_autodetect(_MEMILIO_SETTINGS)
        conan_cmake_install(
            PATH_OR_REFERENCE .
            BUILD missing
            REMOTE conancenter
            SETTINGS ${_MEMILIO_SETTINGS})
    endif()

    list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_BINARY_DIR})
    list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})
endif()

find_package(Boost ${_MEMILIO_BOOST_VERSION} REQUIRED CONFIG COMPONENTS filesystem)

find_package(Eigen3 ${_MEMILIO_EIGEN_VERSION} REQUIRED CONFIG)

find_package(spdlog ${_MEMILIO_SPDLOG_VERSION} CONFIG REQUIRED)

find_package(jsoncpp ${_MEMILIO_JSONCPP_VERSION} CONFIG)
if (TARGET jsoncpp::jsoncpp)
    message(STATUS "MEMILIO WORKAROUND: Creating ALIAS Target JsonCpp::JsonCpp. Target names in the Conan package are lower-case.")
    add_library(JsonCpp::JsonCpp ALIAS jsoncpp::jsoncpp)
endif()
if (TARGET JsonCpp::JsonCpp)
    set(MEMILIO_HAS_JSONCPP ON)
else()
    message(WARNING "JsonCpp was not found. Memilio will be built without some IO features. 
        Set CMake variable MEMILIO_USE_BUNDLED_JSONCPP to ON or install JsonCpp and set the jsoncpp_DIR cmake variable 
        to the directory containing the jsoncppConfig.cmake file to build with JsonCpp.")
endif()

find_package(HDF5 ${_MEMILIO_HDF5_VERSION} CONFIG COMPONENTS C)
if (HDF5_FOUND AND NOT TARGET HDF5::C)
    message(STATUS "MEMILIO WORKAROUND: Creating ALIAS Target HDF5::C. Some HDF5 distributions don't use namespace targets.")
    add_library(HDF5::C INTERFACE)
    target_link_libraries(HDF5::C INTERFACE ${HDF5_C_LIBRARIES})
    target_include_directories(HDF5::C INTERFACE ${HDF5_INCLUDE_DIRS})
endif()
if (TARGET HDF5::C)
    set(MEMILIO_HAS_HDF5 ON)
else()
    message(WARNING "HDF5 was not found. Memilio will be built without some IO features. Install HDF5 Libraries and set the HDF5_DIR cmake variable to the directory containing the hdf5-config.cmake file to build with HDF5.")
endif()