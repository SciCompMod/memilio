cmake_minimum_required(VERSION 3.13)
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

# Function to add a pybind11 module
function(add_pymio_module target_name)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs LINKED_LIBRARIES SOURCES)
    cmake_parse_arguments(PYBIND11_MODULE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})


    pybind11_add_module(${target_name} MODULE ${PYBIND11_MODULE_SOURCES})
    target_link_libraries(${target_name} PRIVATE ${PYBIND11_MODULE_LINKED_LIBRARIES})
    target_include_directories(${target_name} PRIVATE memilio/simulation/bindings)
    install(TARGETS ${target_name} LIBRARY DESTINATION memilio/simulation)
endfunction()


# build python extensions
add_pymio_module(_simulation_abm
    LINKED_LIBRARIES memilio abm
    SOURCES memilio/simulation/bindings/models/abm.cpp
)

add_pymio_module(_simulation
    LINKED_LIBRARIES memilio
    SOURCES memilio/simulation/bindings/simulation.cpp
    memilio/simulation/bindings/epidemiology/damping_sampling.cpp
    memilio/simulation/bindings/epidemiology/uncertain_matrix.cpp
    memilio/simulation/bindings/mobility/metapopulation_mobility_instant.cpp
    memilio/simulation/bindings/utils/date.cpp
    memilio/simulation/bindings/utils/logging.cpp
    memilio/simulation/bindings/utils/time_series.cpp
    memilio/simulation/bindings/utils/parameter_distributions.cpp
    memilio/simulation/bindings/utils/uncertain_value.cpp
)

add_pymio_module(_simulation_osir
    LINKED_LIBRARIES memilio ode_sir
    SOURCES memilio/simulation/bindings/models/osir.cpp
)

add_pymio_module(_simulation_oseir
    LINKED_LIBRARIES memilio ode_seir
    SOURCES memilio/simulation/bindings/models/oseir.cpp
)

add_pymio_module(_simulation_osecir
    LINKED_LIBRARIES memilio ode_secir
    SOURCES memilio/simulation/bindings/models/osecir.cpp
)

add_pymio_module(_simulation_osecirvvs
    LINKED_LIBRARIES memilio ode_secirvvs
    SOURCES memilio/simulation/bindings/models/osecirvvs.cpp
)
