include(FetchContent)
FetchContent_Declare(benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.6.1)

set(BENCHMARK_ENABLE_TESTING OFF)
set(BENCHMARK_ENABLE_INSTALL OFF)
set(BENCHMARK_DOWNLOAD_DEPENDENCIES OFF)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)

FetchContent_MakeAvailable(benchmark)

add_executable(integrator_step_benchmark integrator_step.cpp)
target_link_libraries(integrator_step_benchmark PRIVATE memilio ode_secir benchmark::benchmark)

add_executable(flow_simulation_ode_secirvvs_benchmark flow_simulation_ode_secirvvs.cpp)
target_link_libraries(flow_simulation_ode_secirvvs_benchmark PRIVATE memilio ode_secirvvs benchmark::benchmark)

add_executable(flow_simulation_ode_seir_benchmark flow_simulation_ode_seir.cpp)
target_link_libraries(flow_simulation_ode_seir_benchmark PRIVATE memilio ode_seir benchmark::benchmark)

add_executable(simulation_benchmark simulation.cpp)
target_link_libraries(simulation_benchmark PRIVATE memilio ode_secir benchmark::benchmark)
add_executable(graph_simulation_benchmark graph_simulation.cpp)
target_link_libraries(graph_simulation_benchmark PRIVATE memilio ode_secirvvs benchmark::benchmark)

add_executable(abm_benchmark abm.cpp)
target_link_libraries(abm_benchmark PRIVATE abm benchmark::benchmark)
