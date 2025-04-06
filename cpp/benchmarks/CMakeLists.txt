# Attempt to find pre-installed Google Benchmark
set(MEMILIO_BENCHMARK_VERSION "v1.6.1" CACHE STRING "GoogleTest version to use")

# Configure Google Benchmark build options BEFORE MakeAvailable
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE) # Don't build benchmark's tests
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE) # Don't install benchmark
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE) # Don't require gtest

FetchContent_Declare(google_benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        ${MEMILIO_BENCHMARK_VERSION}
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(google_benchmark)

if(TARGET benchmark::benchmark)
    set(MEMILIO_BENCHMARK_LIBRARY benchmark::benchmark)
else()
    message(FATAL_ERROR "FetchContent for Google Benchmark did not provide target 'benchmark::benchmark'.")
endif() 



# Helper function to add benchmark executables
function(add_memilio_benchmark NAME SOURCES)
    set(MEMILIO_DEPS ${ARGN})

    add_executable(${NAME} ${SOURCES})
    target_link_libraries(${NAME} PRIVATE
        ${MEMILIO_DEPS}
        ${MEMILIO_BENCHMARK_LIBRARY}
    )
    target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}) # For bench_config_dir.h
endfunction()

file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/configs" MEMILIO_BENCHMARK_CONFIG_DIR)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/bench_config_dir.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/bench_config_dir.h"
)

# Define benchmark executables
add_memilio_benchmark(memilio_integrator_step_benchmark
    integrator_step.cpp
    ode_secir
)

add_memilio_benchmark(memilio_flow_simulation_ode_secirvvs_benchmark
    flow_simulation_ode_secirvvs.cpp
    ode_secirvvs
)

add_memilio_benchmark(memilio_flow_simulation_ode_seir_benchmark
    flow_simulation_ode_seir.cpp
    ode_seir
)

add_memilio_benchmark(memilio_simulation_benchmark
    simulation.cpp
    ode_secir
)

add_memilio_benchmark(memilio_graph_simulation_benchmark
    graph_simulation.cpp
    ode_secirvvs
)

add_memilio_benchmark(memilio_abm_benchmark
    abm.cpp
    abm
)