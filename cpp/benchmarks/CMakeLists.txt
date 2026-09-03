set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark testing" FORCE)

# set(BENCHMARK_ENABLE_EXCEPTIONS OFF CACHE BOOL "Disable benchmark exceptions" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Don't install benchmark" FORCE)
set(BENCHMARK_DOWNLOAD_DEPENDENCIES OFF CACHE BOOL "Don't download dependencies" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Disable Google Test in benchmark" FORCE)

message(STATUS "Downloading GoogleBenchmark library")
include(FetchContent)
FetchContent_Declare(benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.6.1)
FetchContent_MakeAvailable(benchmark)

set_target_properties(benchmark PROPERTIES FOLDER "Extern")

add_executable(integrator_step_benchmark integrator_step.cpp)
target_link_libraries(integrator_step_benchmark PRIVATE memilio ode_secir benchmark::benchmark)

add_executable(flow_simulation_ode_secirvvs_benchmark flow_simulation_ode_secirvvs.cpp)
target_link_libraries(flow_simulation_ode_secirvvs_benchmark PRIVATE memilio ode_secirvvs benchmark::benchmark)

add_executable(flow_simulation_ode_seir_benchmark flow_simulation_ode_seir.cpp)
target_link_libraries(flow_simulation_ode_seir_benchmark PRIVATE memilio ode_seir benchmark::benchmark)

add_executable(ode_seir_mobility_outcomes ode_seir_mobility_outcomes.cpp)
target_link_libraries(ode_seir_mobility_outcomes PRIVATE memilio ode_seir_metapop)

add_executable(ode_seir_benchmark_stage_aligned
    ode_seir_benchmark_stage_aligned.cpp
    ode_seir_benchmark_stage_aligned.h)
target_link_libraries(ode_seir_benchmark_stage_aligned PRIVATE memilio benchmark::benchmark)
if(MEMILIO_ENABLE_OPENMP)
    target_link_libraries(ode_seir_benchmark_stage_aligned PRIVATE OpenMP::OpenMP_CXX)
endif()
if(MEMILIO_NATIVE_BENCHMARKS)
    target_compile_options(ode_seir_benchmark_stage_aligned PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-march=native>)
endif()

add_executable(ode_seir_metapop_benchmark_cpu
    ode_seir_metapop_benchmark_cpu.cpp
    ode_seir_metapop_benchmark.h)
target_link_libraries(ode_seir_metapop_benchmark_cpu PRIVATE Eigen3::Eigen benchmark::benchmark)
target_compile_definitions(ode_seir_metapop_benchmark_cpu PRIVATE EIGEN_DONT_PARALLELIZE)
if(MEMILIO_ENABLE_OPENMP)
    target_link_libraries(ode_seir_metapop_benchmark_cpu PRIVATE OpenMP::OpenMP_CXX)
endif()
if(MEMILIO_NATIVE_BENCHMARKS)
    target_compile_options(ode_seir_metapop_benchmark_cpu PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang>:-march=native>)
endif()

if(MEMILIO_BUILD_CUDA_BENCHMARKS)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    add_executable(ode_seir_benchmark_gpu
        ode_seir_benchmark_gpu.cpp
        ode_seir_benchmark_gpu_kernels.cu
        ode_seir_benchmark_stage_aligned.h)
    target_link_libraries(ode_seir_benchmark_gpu PRIVATE benchmark::benchmark CUDA::cudart)
    set_target_properties(ode_seir_benchmark_gpu PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION OFF)

    add_executable(ode_seir_metapop_benchmark_gpu
        ode_seir_metapop_benchmark_gpu.cu
        ode_seir_metapop_benchmark.h)
    target_link_libraries(ode_seir_metapop_benchmark_gpu PRIVATE
        Eigen3::Eigen benchmark::benchmark CUDA::cudart CUDA::cublas)
    target_compile_definitions(ode_seir_metapop_benchmark_gpu PRIVATE EIGEN_DONT_PARALLELIZE EIGEN_NO_CUDA)
    set_target_properties(ode_seir_metapop_benchmark_gpu PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION OFF)
    if(MEMILIO_NATIVE_BENCHMARKS)
        target_compile_options(ode_seir_metapop_benchmark_gpu PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-march=native>)
    endif()
endif()

add_executable(simulation_benchmark simulation.cpp)
target_link_libraries(simulation_benchmark PRIVATE memilio ode_secir benchmark::benchmark)
add_executable(graph_simulation_benchmark graph_simulation.cpp)
target_link_libraries(graph_simulation_benchmark PRIVATE memilio ode_secirvvs benchmark::benchmark)

add_executable(secirvvs_advance_benchmark secirvvs_advance.cpp)
target_link_libraries(secirvvs_advance_benchmark PRIVATE memilio ode_secirvvs benchmark::benchmark)

add_executable(abm_benchmark abm.cpp)
target_link_libraries(abm_benchmark PRIVATE abm benchmark::benchmark)
