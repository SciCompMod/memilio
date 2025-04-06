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
        memilio
        ${MEMILIO_DEPS}
        ${MEMILIO_BENCHMARK_LIBRARY}
    )
    memilio_configure_path_header(${NAME} "bench_config_dir.h.in" "bench_config_dir.h" "MEMILIO_BENCHMARK_CONFIG_DIR" "${CMAKE_CURRENT_SOURCE_DIR}/configs")
endfunction()


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