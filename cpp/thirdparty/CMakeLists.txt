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
    message(STATUS "Memilio: Installing required dependencies with Conan.")
    
    include(memilio-conan)

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

    memilio_conan_install(
        REQUIRES ${_MEMILIO_REQUIRES}
        OPTIONS
            boost:without_atomic=True
            boost:without_chrono=True
            boost:without_container=True
            boost:without_context=True
            boost:without_contract=True
            boost:without_coroutine=True
            boost:without_date_time=True
            boost:without_exception=True
            boost:without_fiber=True
            boost:without_filesystem=False
            boost:without_graph=True
            boost:without_graph_parallel=True
            boost:without_iostreams=True
            boost:without_json=True
            boost:without_locale=True
            boost:without_log=True
            boost:without_math=True
            boost:without_mpi=True
            boost:without_nowide=True
            boost:without_program_options=True
            boost:without_python=True
            boost:without_random=True
            boost:without_regex=True
            boost:without_serialization=True
            boost:without_stacktrace=True
            boost:without_system=False
            boost:without_test=True
            boost:without_thread=True
            boost:without_timer=True
            boost:without_type_erasure=True
            boost:without_wave=True)

    # force find_package to use the Conan packages instead of other distributions on the system
    if (NOT MEMILIO_USE_SYSTEM_BOOST)
        set(_MEMILIO_BOOST_IMPORT CONFIG)
        set(Boost_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
    endif()
    if (NOT MEMILIO_USE_SYSTEM_EIGEN)
        set(_MEMILIO_EIGEN_IMPORT CONFIG)
        set(Eigen3_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
    endif()
    if (NOT MEMILIO_USE_SYSTEM_SPDLOG)
        set(_MEMILIO_SPDLOG_IMPORT CONFIG)
        set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
    endif()
    if (NOT MEMILIO_USE_SYSTEM_JSONCPP)
        set(_MEMILIO_JSONCPP_IMPORT CONFIG)
        set(jsoncpp_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
    endif()
    if (NOT MEMILIO_USE_SYSTEM_HDF5)
        set(_MEMILIO_HDF5_IMPORT CONFIG)
        set(HDF5_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
    endif()
endif()

find_package(Boost REQUIRED ${_MEMILIO_BOOST_IMPORT} COMPONENTS filesystem)

find_package(Eigen3 REQUIRED ${_MEMILIO_EIGEN_IMPORT})

find_package(spdlog REQUIRED ${_MEMILIO_SPDLOG_IMPORT})

find_package(jsoncpp ${_MEMILIO_JSONCPP_IMPORT})
memilio_target_alternatives(TARGET JsonCpp::JsonCpp ALTERNATIVES jsoncpp::jsoncpp)
if (TARGET JsonCpp::JsonCpp)
    set(MEMILIO_HAS_JSONCPP ON)
else()
    message(WARNING "JsonCpp was not found. Memilio will be built without some IO features. 
        Set CMake variable MEMILIO_USE_BUNDLED_JSONCPP to ON or install JsonCpp and set the jsoncpp_DIR cmake variable 
        to the directory containing the jsoncppConfig.cmake file to build with JsonCpp.")
endif()

find_package(HDF5 ${_MEMILIO_HDF5_IMPORT} COMPONENTS C)
if (HDF5_FOUND AND NOT TARGET HDF5::C)
    message(STATUS "MEMILIO WORKAROUND: Creating ALIAS Target HDF5::C. Some HDF5 distributions don't use namespace targets.")
    add_library(HDF5::C INTERFACE IMPORTED)
    target_link_libraries(HDF5::C INTERFACE ${HDF5_C_LIBRARIES})
    target_include_directories(HDF5::C INTERFACE ${HDF5_INCLUDE_DIRS})
endif()
if (TARGET HDF5::C)
    set(MEMILIO_HAS_HDF5 ON)
else()
    message(WARNING "HDF5 was not found. Memilio will be built without some IO features. Install HDF5 Libraries and set the HDF5_DIR cmake variable to the directory containing the hdf5-config.cmake file to build with HDF5.")
endif()
