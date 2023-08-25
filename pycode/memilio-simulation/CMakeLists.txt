cmake_minimum_required(VERSION 3.4)
project(memilio-python)

option(MEMILIO_USE_BUNDLED_PYBIND11 "Use pybind11 bundled with this library." ON)
mark_as_advanced(MEMILIO_USE_BUNDLED_PYBIND11)

include(GNUInstallDirs) # set to gnu folders. No cache variable so this is not global

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_RPATH "${CMAKE_BINARY_DIR}/lib" "${CMAKE_BINARY_DIR}/bin")

if(MEMILIO_USE_BUNDLED_PYBIND11)
    # Fetch pybind11
    include(FetchContent)
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG v2.11.1
    )
    FetchContent_GetProperties(pybind11)

    if(NOT pybind11_POPULATED)
        FetchContent_Populate(pybind11)
        add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif()
else()
    find_package(pybind11 REQUIRED)
endif()

# add in C++ library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../cpp ${CMAKE_CURRENT_BINARY_DIR}/cpp EXCLUDE_FROM_ALL)

# build python extensions
pybind11_add_module(_simulation_secir MODULE
    memilio/simulation/secir.cpp
)
target_link_libraries(_simulation_secir PRIVATE memilio ode_secir)
target_include_directories(_simulation_secir PRIVATE memilio/simulation)
install(TARGETS _simulation_secir LIBRARY DESTINATION memilio)

pybind11_add_module(_simulation_abm MODULE
    memilio/simulation/abm.cpp
)
target_link_libraries(_simulation_abm PRIVATE memilio abm)
target_include_directories(_simulation_abm PRIVATE memilio/simulation)
install(TARGETS _simulation_abm LIBRARY DESTINATION memilio)

pybind11_add_module(_simulation MODULE
    memilio/simulation/simulation.cpp
    memilio/simulation/epidemiology/damping_sampling.cpp
    memilio/simulation/epidemiology/uncertain_matrix.cpp
    memilio/simulation/mobility/metapopulation_mobility_instant.cpp
    memilio/simulation/utils/date.cpp
    memilio/simulation/utils/logging.cpp
    memilio/simulation/utils/time_series.cpp
    memilio/simulation/utils/parameter_distributions.cpp
    memilio/simulation/utils/uncertain_value.cpp
)
target_link_libraries(_simulation PRIVATE memilio)
target_include_directories(_simulation PRIVATE memilio/simulation)
install(TARGETS _simulation LIBRARY DESTINATION memilio)

pybind11_add_module(_simulation_oseir MODULE
    memilio/simulation/oseir.cpp
)
target_link_libraries(_simulation_oseir PRIVATE memilio ode_seir)
target_include_directories(_simulation_oseir PRIVATE memilio/simulation)
install(TARGETS _simulation_oseir LIBRARY DESTINATION memilio)
