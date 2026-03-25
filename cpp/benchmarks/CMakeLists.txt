set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark testing" FORCE)

# set(BENCHMARK_ENABLE_EXCEPTIONS OFF CACHE BOOL "Disable benchmark exceptions" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Don't install benchmark" FORCE)
set(BENCHMARK_DOWNLOAD_DEPENDENCIES OFF CACHE BOOL "Don't download dependencies" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Disable Google Test in benchmark" FORCE)

include(FetchContent)
FetchContent_Declare(benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.6.1)
FetchContent_GetProperties(benchmark)

if(NOT benchmark_POPULATED)
    FetchContent_Populate(benchmark)
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE BOOL "")
    add_subdirectory(${benchmark_SOURCE_DIR} ${benchmark_BINARY_DIR} EXCLUDE_FROM_ALL)
    unset(CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
endif()

set_target_properties(benchmark PROPERTIES FOLDER "Extern")

add_executable(integrator_step_benchmark integrator_step.cpp)
target_link_libraries(integrator_step_benchmark PRIVATE memilio ode_secir benchmark::benchmark)

add_executable(ode_seir_benchmark ode_seir_benchmark.cpp)
target_link_libraries(ode_seir_benchmark PRIVATE memilio ode_seir benchmark::benchmark)

add_executable(ode_seir_benchmark_w_sim ode_seir_benchmark_w_sim.cpp)
target_link_libraries(ode_seir_benchmark_w_sim PRIVATE memilio ode_seir benchmark::benchmark)

option(WITH_LIKWID "Enable LIKWID Marker API instrumentation" OFF)

if(WITH_LIKWID)
    find_path(LIKWID_INCLUDE_DIR likwid.h
        HINTS /usr/include /usr/local/include)
    find_library(LIKWID_LIB likwid
        HINTS /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu)

    if(LIKWID_INCLUDE_DIR AND LIKWID_LIB)
        message(STATUS "LIKWID found: ${LIKWID_LIB} – enabling Marker API")
        target_compile_definitions(ode_seir_benchmark_w_sim PRIVATE LIKWID_PERFMON)
        target_include_directories(ode_seir_benchmark_w_sim PRIVATE "${LIKWID_INCLUDE_DIR}")
        target_link_libraries(ode_seir_benchmark_w_sim PRIVATE "${LIKWID_LIB}")
    else()
        message(WARNING "WITH_LIKWID=ON but likwid.h / liblikwid not found – disabling")
    endif()
endif()

find_package(OpenMP QUIET)

if(OpenMP_CXX_FOUND)
    target_link_libraries(ode_seir_benchmark_w_sim PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP found – bench_stage_aligned_rk4_omp will be compiled")
else()
    message(STATUS "OpenMP NOT found – bench_stage_aligned_rk4_omp will be skipped")
endif()

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

add_executable(ode_secirts_benchmark_stage_aligned ode_secirts_benchmark_stage_aligned.cpp)
target_link_libraries(ode_secirts_benchmark_stage_aligned PRIVATE memilio ode_secirts benchmark::benchmark)

# Optional: GPU-accelerated stage-aligned commuter benchmark (requires CUDA)
# Enable only when a CUDA compiler is available.  To configure explicitly:
# cmake -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc ..
include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    get_filename_component(_cuda_bin_dir "${CMAKE_CUDA_COMPILER}" DIRECTORY)
    get_filename_component(_cuda_root "${_cuda_bin_dir}" DIRECTORY)
    find_library(_cuda_rt_lib cudart
        HINTS "${_cuda_root}/lib64" "${_cuda_root}/lib"
        "${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES}"
        NO_DEFAULT_PATH)

    if(NOT _cuda_rt_lib)
        find_library(_cuda_rt_lib cudart)
    endif()

    add_library(seir_gpu_kernels OBJECT
        ode_seir_benchmark_gpu_kernels.cu)
    set_target_properties(seir_gpu_kernels PROPERTIES
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        CUDA_SEPARABLE_COMPILATION OFF
        POSITION_INDEPENDENT_CODE ON
    )
    target_compile_options(seir_gpu_kernels PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math;-O3;--expt-relaxed-constexpr>
    )

    add_executable(ode_seir_benchmark_gpu
        ode_seir_benchmark_gpu.cpp
        $<TARGET_OBJECTS:seir_gpu_kernels>)
    target_link_libraries(ode_seir_benchmark_gpu
        PRIVATE memilio ode_seir benchmark::benchmark "${_cuda_rt_lib}")
    target_include_directories(ode_seir_benchmark_gpu PRIVATE
        "${_cuda_root}/include")
    target_compile_options(ode_seir_benchmark_gpu PRIVATE -O3)

    message(STATUS "CUDA found (${CMAKE_CUDA_COMPILER}) – "
        "building ode_seir_benchmark_gpu")
    message(STATUS "cudart: ${_cuda_rt_lib}")
else()
    message(STATUS "No CUDA compiler found – "
        "skipping ode_seir_benchmark_gpu")
endif()
