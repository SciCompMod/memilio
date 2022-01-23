include(memilio-conan)

message(STATUS "Memilio: Installing required dependencies with Conan.")

memilio_conan_install(
    REQUIRES 
        boost/1.75.0
        eigen/3.3.9
        spdlog/1.5.0
        jsoncpp/1.9.5
        hdf5/1.12.0
    OPTIONS
        *:shared=False
        *:fPIC=True
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
        boost:without_wave=True
        hdf5:enable_cxx=False
        hdf5:hl=False)

set(Boost_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
set(Eigen3_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
set(spdlog_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
set(jsoncpp_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)
set(HDF5_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "" FORCE)

find_package(Boost REQUIRED CONFIG COMPONENTS filesystem)
find_package(Eigen3 REQUIRED CONFIG)
find_package(spdlog REQUIRED CONFIG)
find_package(jsoncpp REQUIRED CONFIG)
set(MEMILIO_HAS_JSONCPP ON)
find_package(HDF5 REQUIRED CONFIG COMPONENTS C)
set(MEMILIO_HAS_HDF5 ON)
