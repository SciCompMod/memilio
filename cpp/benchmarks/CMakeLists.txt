set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark testing" FORCE)

# set(BENCHMARK_ENABLE_EXCEPTIONS OFF CACHE BOOL "Disable benchmark exceptions" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Don't install benchmark" FORCE)
set(BENCHMARK_DOWNLOAD_DEPENDENCIES OFF CACHE BOOL "Don't download dependencies" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Disable Google Test in benchmark" FORCE)


include(FetchContent)
FetchContent_Declare(benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.6.1)
FetchContent_MakeAvailable(benchmark)


set_target_properties(benchmark PROPERTIES FOLDER "Extern")

add_executable(simulation_benchmark simulation.cpp)
target_link_libraries(simulation_benchmark PRIVATE memilio ode_secir benchmark::benchmark)

add_executable(integrator_step_benchmark integrator_step.cpp)
target_link_libraries(integrator_step_benchmark PRIVATE memilio ode_secir benchmark::benchmark)
