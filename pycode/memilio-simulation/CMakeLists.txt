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

# Function to add a pybind11 module
function(add_pymio_module target_name)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs LINKED_LIBRARIES SOURCES)
    cmake_parse_arguments(PYBIND11_MODULE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})


    pybind11_add_module(${target_name} MODULE ${PYBIND11_MODULE_SOURCES})
    target_link_libraries(${target_name} PRIVATE ${PYBIND11_MODULE_LINKED_LIBRARIES})
    target_include_directories(${target_name} PRIVATE memilio/simulation)
    install(TARGETS ${target_name} LIBRARY DESTINATION memilio)
endfunction()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(memilio/simulation)

# # build python extensions
add_library(osir_obj OBJECT memilio/simulation/models/osir.cpp)
target_link_libraries(osir_obj PRIVATE memilio ode_sir pybind11::module)

add_library(oseir_obj OBJECT memilio/simulation/models/oseir.cpp)
target_link_libraries(oseir_obj PRIVATE memilio ode_seir pybind11::module)

add_library(osecir_obj OBJECT memilio/simulation/models/osecir.cpp)
target_link_libraries(osecir_obj PRIVATE memilio ode_secir pybind11::module)

add_library(osecirvvs_obj OBJECT memilio/simulation/models/osecirvvs.cpp)
target_link_libraries(osecirvvs_obj PRIVATE memilio ode_secirvvs pybind11::module)

add_library(abm_obj OBJECT memilio/simulation/models/abm.cpp)
target_link_libraries(abm_obj PRIVATE memilio abm pybind11::module)

add_library(simulation_obj OBJECT 
    memilio/simulation/simulation.cpp
    memilio/simulation/epidemiology/damping_sampling.cpp
    memilio/simulation/epidemiology/uncertain_matrix.cpp
    memilio/simulation/mobility/metapopulation_mobility_instant.cpp
    memilio/simulation/utils/date.cpp
    memilio/simulation/utils/logging.cpp
    memilio/simulation/utils/time_series.cpp
    memilio/simulation/utils/parameter_distributions.cpp
    memilio/simulation/utils/uncertain_value.cpp)
target_link_libraries(simulation_obj PRIVATE memilio pybind11::module)


# add_pymio_module(simulation
#     LINKED_LIBRARIES osir_obj
#     oseir_obj
#     osecir_obj
#     osecirvvs_obj
#     abm_obj
#     simulation_obj
#     SOURCES memilio/simulation/pymio.cpp
# )

add_pymio_module(simulation
    LINKED_LIBRARIES memilio ode_sir ode_seir ode_secir ode_secirvvs abm
    SOURCES memilio/simulation/pymio.cpp
    $<TARGET_OBJECTS:osir_obj>
    $<TARGET_OBJECTS:oseir_obj>
    $<TARGET_OBJECTS:osecir_obj>
    $<TARGET_OBJECTS:osecirvvs_obj>
    $<TARGET_OBJECTS:abm_obj>
    $<TARGET_OBJECTS:simulation_obj>
)