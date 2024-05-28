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

# Find the Python interpreter
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Specify the output directory for the stubs
set(STUBS_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/memilio-stubs/simulation/")

# build python extensions
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
add_custom_command(
    TARGET _simulation
    POST_BUILD
    COMMAND ${Python3_EXECUTABLE} -m pybind11_stubgen -o ${STUBS_OUTPUT_DIR} memilio._simulation
    COMMAND ${CMAKE_COMMAND} -E copy ${STUBS_OUTPUT_DIR}/memilio/_simulation.pyi ${STUBS_OUTPUT_DIR}/__init__.pyi
    COMMENT "Generating Python stubs for _simulation"
)
install(TARGETS _simulation LIBRARY DESTINATION memilio)

pybind11_add_module(_simulation_osir MODULE
    memilio/simulation/osir.cpp
)
target_link_libraries(_simulation_osir PRIVATE memilio ode_sir)
target_include_directories(_simulation_osir PRIVATE memilio/simulation)
add_custom_command(
    TARGET _simulation_osir
    POST_BUILD
    COMMAND ${Python3_EXECUTABLE} -m pybind11_stubgen -o ${STUBS_OUTPUT_DIR} memilio._simulation_osir
    COMMAND ${CMAKE_COMMAND} -E copy ${STUBS_OUTPUT_DIR}/memilio/_simulation_osir.pyi ${STUBS_OUTPUT_DIR}/osir.pyi
    COMMENT "Generating Python stubs for _simulation_osir"
)
install(TARGETS _simulation_osir LIBRARY DESTINATION memilio)

pybind11_add_module(_simulation_oseir MODULE
    memilio/simulation/oseir.cpp
)
target_link_libraries(_simulation_oseir PRIVATE memilio ode_seir)
target_include_directories(_simulation_oseir PRIVATE memilio/simulation)
# Custom command to generate stubs
add_custom_command(
    TARGET _simulation_oseir
    POST_BUILD
    COMMAND ${Python3_EXECUTABLE} -m pybind11_stubgen -o ${STUBS_OUTPUT_DIR} memilio._simulation_oseir
    COMMAND ${CMAKE_COMMAND} -E copy ${STUBS_OUTPUT_DIR}/memilio/_simulation_oseir.pyi ${STUBS_OUTPUT_DIR}/oseir.pyi
    COMMENT "Generating Python stubs for _simulation_oseir"
)
install(TARGETS _simulation_oseir LIBRARY DESTINATION memilio)

pybind11_add_module(_simulation_secir MODULE
    memilio/simulation/secir.cpp
)
target_link_libraries(_simulation_secir PRIVATE memilio ode_secir)
target_include_directories(_simulation_secir PRIVATE memilio/simulation)
# Custom command to generate stubs
add_custom_command(
    TARGET _simulation_secir
    POST_BUILD
    COMMAND ${Python3_EXECUTABLE} -m pybind11_stubgen -o ${STUBS_OUTPUT_DIR} memilio._simulation_secir
    COMMAND ${CMAKE_COMMAND} -E copy ${STUBS_OUTPUT_DIR}/memilio/_simulation_secir.pyi ${STUBS_OUTPUT_DIR}/secir.pyi
    COMMENT "Generating Python stubs for _simulation_secir"
)
install(TARGETS _simulation_secir LIBRARY DESTINATION memilio)

pybind11_add_module(_simulation_osecirvvs MODULE
    memilio/simulation/osecirvvs.cpp
)
target_link_libraries(_simulation_osecirvvs PRIVATE memilio ode_secirvvs)
target_include_directories(_simulation_osecirvvs PRIVATE memilio/simulation)
add_custom_command(
    TARGET _simulation_osecirvvs
    POST_BUILD
    COMMAND ${Python3_EXECUTABLE} -m pybind11_stubgen -o ${STUBS_OUTPUT_DIR} memilio._simulation_osecirvvs
    COMMAND ${CMAKE_COMMAND} -E copy ${STUBS_OUTPUT_DIR}/memilio/_simulation_osecirvvs.pyi ${STUBS_OUTPUT_DIR}/osecirvvs.pyi
    COMMENT "Generating Python stubs for _simulation_osecirvvs"
)
install(TARGETS _simulation_osecirvvs LIBRARY DESTINATION memilio)