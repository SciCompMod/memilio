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
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS TRUE)
FetchContent_MakeAvailable(google_benchmark)
unset(CMAKE_SUPPRESS_DEVELOPER_WARNINGS)


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
endfunction()


# Define benchmark executables
add_memilio_benchmark(memilio_abm_benchmark
    abm.cpp
    abm
)

if(MEMILIO_HAS_JSONCPP)
    add_memilio_benchmark(memilio_graph_simulation_benchmark
        graph_simulation.cpp
        ode_secirvvs
    )

    add_memilio_benchmark(memilio_flow_simulation_ode_seir_benchmark
        flow_simulation_ode_seir.cpp
        ode_seir
    )


    add_memilio_benchmark(memilio_flow_simulation_ode_secirvvs_benchmark
        flow_simulation_ode_secirvvs.cpp
        ode_secirvvs
    )

    add_memilio_benchmark(integrator_secir_benchmark
        integrator_secir_sim.cpp
        ode_secir
    )

    add_memilio_benchmark(memilio_integrator_step_benchmark
        integrator_step.cpp
        ode_secir
    )
endif()