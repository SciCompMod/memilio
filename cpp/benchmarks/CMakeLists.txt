include(FetchContent)
FetchContent_Declare(benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.6.1)

set(BENCHMARK_ENABLE_TESTING OFF)
set(BENCHMARK_ENABLE_INSTALL OFF)
set(BENCHMARK_DOWNLOAD_DEPENDENCIES OFF)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)

FetchContent_MakeAvailable(benchmark)

add_executable(simulation_benchmark simulation.cpp)
target_link_libraries(simulation_benchmark PRIVATE memilio ode_secir benchmark::benchmark)

add_executable(integrator_step_benchmark integrator_step.cpp)
target_link_libraries(integrator_step_benchmark PRIVATE memilio ode_secir benchmark::benchmark)
