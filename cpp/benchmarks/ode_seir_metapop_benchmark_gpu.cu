/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#include "benchmark/benchmark.h"
#include "ode_seir_metapop_benchmark.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mio::benchmark_mio
{

void check_cuda(cudaError_t status, const char* operation)
{
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
    }
}

void check_cublas(cublasStatus_t status, const char* operation)
{
    if (status != CUBLAS_STATUS_SUCCESS) {
        throw std::runtime_error(std::string(operation) + " failed with cuBLAS status " +
                                 std::to_string(static_cast<int>(status)) + ".");
    }
}

__constant__ double device_infection_coefficients[maximum_groups * maximum_groups];
__constant__ double device_rate_exposed[maximum_groups];
__constant__ double device_rate_infected[maximum_groups];

__global__ void normalize_present_share(double* infectious_share, const double* inverse_present_population,
                                        size_t values)
{
    const size_t index = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index < values) {
        infectious_share[index] *= inverse_present_population[index];
    }
}

template <int Stage>
__device__ __forceinline__ void update_rk_value_device(double base, double derivative, double dt, double& accumulated,
                                                       double& stage_value, double& final_value)
{
    if constexpr (Stage == 0) {
        accumulated = fma(dt / 6.0, derivative, base);
        stage_value = fma(dt / 2.0, derivative, base);
    }
    else if constexpr (Stage == 1) {
        accumulated = fma(dt / 3.0, derivative, accumulated);
        stage_value = fma(dt / 2.0, derivative, base);
    }
    else if constexpr (Stage == 2) {
        accumulated = fma(dt / 3.0, derivative, accumulated);
        stage_value = fma(dt, derivative, base);
    }
    else {
        final_value = fma(dt / 6.0, derivative, accumulated);
    }
}

template <int G, int Stage>
__global__ void rhs_and_rk_kernel(double* __restrict__ state, double* __restrict__ stage_state,
                                  double* __restrict__ result, const double* __restrict__ mobile_share,
                                  const double* __restrict__ inverse_population, int patches, double dt)
{
    const int patch = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (patch >= patches) {
        return;
    }
    const double* current = Stage == 0 ? state : stage_state;
    double infectious_share[G];
#pragma unroll
    for (int source = 0; source < G; ++source) {
        const size_t infected = static_cast<size_t>(3 * source + 2) * patches + patch;
        infectious_share[source] =
            current[infected] * inverse_population[static_cast<size_t>(source) * patches + patch];
    }
#pragma unroll
    for (int target = 0; target < G; ++target) {
        double pressure = 0.0;
#pragma unroll
        for (int source = 0; source < G; ++source) {
            pressure =
                fma(device_infection_coefficients[target * G + source],
                    infectious_share[source] + mobile_share[static_cast<size_t>(source) * patches + patch], pressure);
        }
        const size_t susceptible = static_cast<size_t>(3 * target) * patches + patch;
        const size_t exposed     = susceptible + patches;
        const size_t infected    = exposed + patches;
        const double flow_se     = pressure * current[susceptible];
        const double flow_ei     = device_rate_exposed[target] * current[exposed];
        const double flow_ir     = device_rate_infected[target] * current[infected];
        update_rk_value_device<Stage>(state[susceptible], -flow_se, dt, result[susceptible], stage_state[susceptible],
                                      state[susceptible]);
        update_rk_value_device<Stage>(state[exposed], flow_se - flow_ei, dt, result[exposed], stage_state[exposed],
                                      state[exposed]);
        update_rk_value_device<Stage>(state[infected], flow_ei - flow_ir, dt, result[infected], stage_state[infected],
                                      state[infected]);
    }
}

template <int G>
class GpuRunner
{
public:
    explicit GpuRunner(const ImplicitProblem& problem)
        : m_patches(problem.patches)
        , m_state_size(problem.state.size())
    {
        try {
            check_cuda(cudaStreamCreateWithFlags(&m_stream, cudaStreamNonBlocking), "cudaStreamCreateWithFlags");
            check_cublas(cublasCreate(&m_cublas), "cublasCreate");
            check_cublas(cublasSetStream(m_cublas, m_stream), "cublasSetStream");
            allocate(m_state, problem.state.size());
            allocate(m_stage_state, problem.stage_state.size());
            allocate(m_result, problem.result.size());
            allocate(m_present_share, problem.present_share.size());
            allocate(m_mobile_share, problem.mobile_share.size());
            allocate(m_inverse_population, problem.inverse_population.size());
            allocate(m_inverse_present_population, problem.inverse_present_population.size());
            allocate(m_commuting, problem.commuting_row_major.size());
            copy_to_device(m_inverse_population, problem.inverse_population);
            copy_to_device(m_inverse_present_population, problem.inverse_present_population);
            copy_to_device(m_commuting, problem.commuting_row_major);
            check_cuda(cudaMemcpyToSymbol(device_infection_coefficients, problem.infection_coefficients.data(),
                                          problem.infection_coefficients.size() * sizeof(double)),
                       "copy infection coefficients");
            check_cuda(cudaMemcpyToSymbol(device_rate_exposed, problem.rate_exposed.data(),
                                          problem.rate_exposed.size() * sizeof(double)),
                       "copy exposed rates");
            check_cuda(cudaMemcpyToSymbol(device_rate_infected, problem.rate_infected.data(),
                                          problem.rate_infected.size() * sizeof(double)),
                       "copy infected rates");
            check_cuda(cudaMemsetAsync(m_state, 0, m_state_size * sizeof(double), m_stream), "initialize state");
            enqueue_step();
            check_cuda(cudaStreamSynchronize(m_stream), "CUDA warm-up");
            capture_graph();
        }
        catch (...) {
            release();
            throw;
        }
    }

    ~GpuRunner()
    {
        release();
    }

    GpuRunner(const GpuRunner&)            = delete;
    GpuRunner& operator=(const GpuRunner&) = delete;

    void reset(const std::vector<double>& state)
    {
        check_cuda(
            cudaMemcpyAsync(m_state, state.data(), state.size() * sizeof(double), cudaMemcpyHostToDevice, m_stream),
            "upload state");
        check_cuda(cudaStreamSynchronize(m_stream), "finish state upload");
    }

    void run()
    {
        check_cuda(cudaGraphLaunch(m_graph_exec, m_stream), "cudaGraphLaunch");
        check_cuda(cudaStreamSynchronize(m_stream), "finish CUDA integration");
    }

    std::vector<double> download() const
    {
        std::vector<double> result(m_state_size);
        check_cuda(cudaMemcpy(result.data(), m_state, result.size() * sizeof(double), cudaMemcpyDeviceToHost),
                   "download state");
        return result;
    }

    const double* device_state() const
    {
        return m_state;
    }

private:
    template <class T>
    void allocate(T*& pointer, size_t count)
    {
        check_cuda(cudaMalloc(reinterpret_cast<void**>(&pointer), count * sizeof(T)), "cudaMalloc");
    }

    void copy_to_device(double* destination, const std::vector<double>& source)
    {
        check_cuda(cudaMemcpyAsync(destination, source.data(), source.size() * sizeof(double), cudaMemcpyHostToDevice,
                                   m_stream),
                   "upload constant data");
    }

    template <int Stage>
    void enqueue_stage()
    {
        const double* current = Stage == 0 ? m_state : m_stage_state;
        // m_commuting stores row-major H, which cuBLAS views as column-major H^T.
        check_cublas(cublasDgemm(m_cublas, CUBLAS_OP_N, CUBLAS_OP_N, m_patches, G, m_patches, &m_one,
                                 m_commuting, m_patches, current + 2 * m_patches, 3 * m_patches, &m_zero,
                                 m_present_share, m_patches),
                     "cublasDgemm(H^T I)");
        constexpr int point_threads = 256;
        const size_t share_values   = static_cast<size_t>(G) * m_patches;
        const int share_blocks      = static_cast<int>((share_values + point_threads - 1) / point_threads);
        normalize_present_share<<<share_blocks, point_threads, 0, m_stream>>>(
            m_present_share, m_inverse_present_population, share_values);
        check_cuda(cudaGetLastError(), "normalize_present_share");
        check_cublas(cublasDgemm(m_cublas, CUBLAS_OP_T, CUBLAS_OP_N, m_patches, G, m_patches, &m_one,
                                 m_commuting, m_patches, m_present_share, m_patches, &m_zero, m_mobile_share,
                                 m_patches),
                     "cublasDgemm(H U)");
        const int point_blocks = (m_patches + 255) / 256;
        rhs_and_rk_kernel<G, Stage><<<point_blocks, 256, 0, m_stream>>>(
            m_state, m_stage_state, m_result, m_mobile_share, m_inverse_population, m_patches, step_size);
        check_cuda(cudaGetLastError(), "rhs_and_rk_kernel");
    }

    void enqueue_step()
    {
        enqueue_stage<0>();
        enqueue_stage<1>();
        enqueue_stage<2>();
        enqueue_stage<3>();
    }

    void capture_graph()
    {
        check_cuda(cudaStreamBeginCapture(m_stream, cudaStreamCaptureModeGlobal), "cudaStreamBeginCapture");
        for (int step = 0; step < integration_steps; ++step) {
            enqueue_step();
        }
        check_cuda(cudaStreamEndCapture(m_stream, &m_graph), "cudaStreamEndCapture");
        check_cuda(cudaGraphInstantiate(&m_graph_exec, m_graph, nullptr, nullptr, 0), "cudaGraphInstantiate");
        check_cuda(cudaGraphUpload(m_graph_exec, m_stream), "cudaGraphUpload");
        check_cuda(cudaStreamSynchronize(m_stream), "finish CUDA graph upload");
    }

    void release() noexcept
    {
        if (m_graph_exec != nullptr) {
            cudaGraphExecDestroy(m_graph_exec);
        }
        if (m_graph != nullptr) {
            cudaGraphDestroy(m_graph);
        }
        cudaFree(m_commuting);
        cudaFree(m_inverse_present_population);
        cudaFree(m_inverse_population);
        cudaFree(m_mobile_share);
        cudaFree(m_present_share);
        cudaFree(m_result);
        cudaFree(m_stage_state);
        cudaFree(m_state);
        if (m_cublas != nullptr) {
            cublasDestroy(m_cublas);
        }
        if (m_stream != nullptr) {
            cudaStreamDestroy(m_stream);
        }
    }

    int m_patches;
    size_t m_state_size;
    cudaStream_t m_stream                = nullptr;
    cublasHandle_t m_cublas              = nullptr;
    cudaGraph_t m_graph                  = nullptr;
    cudaGraphExec_t m_graph_exec         = nullptr;
    double* m_state                      = nullptr;
    double* m_stage_state                = nullptr;
    double* m_result                     = nullptr;
    double* m_present_share              = nullptr;
    double* m_mobile_share               = nullptr;
    double* m_inverse_population         = nullptr;
    double* m_inverse_present_population = nullptr;
    double* m_commuting                  = nullptr;
    double m_one                         = 1.0;
    double m_zero                        = 0.0;
};

template <int G>
void benchmark_serial_impl(benchmark::State& state)
{
    ImplicitProblem problem(static_cast<int>(state.range(0)), G);
    for (auto _ : state) {
        state.PauseTiming();
        problem.reset_state();
        state.ResumeTiming();
        advance_cpu<G>(problem);
        benchmark::DoNotOptimize(problem.state.data());
        benchmark::ClobberMemory();
    }
    state.counters["patches"]    = problem.patches;
    state.counters["age_groups"] = G;
    state.counters["steps"]      = integration_steps;
}

void benchmark_serial(benchmark::State& state)
{
    switch (static_cast<int>(state.range(1))) {
    case 1:
        benchmark_serial_impl<1>(state);
        break;
    case 3:
        benchmark_serial_impl<3>(state);
        break;
    case 6:
        benchmark_serial_impl<6>(state);
        break;
    case 8:
        benchmark_serial_impl<8>(state);
        break;
    default:
        state.SkipWithError("Unsupported age-group count.");
    }
}

template <int G>
void benchmark_cuda_impl(benchmark::State& state)
{
    try {
        ImplicitProblem problem(static_cast<int>(state.range(0)), G);
        GpuRunner<G> runner(problem);
        for (auto _ : state) {
            state.PauseTiming();
            runner.reset(problem.initial_state);
            state.ResumeTiming();
            runner.run();
            benchmark::DoNotOptimize(runner.device_state());
        }
        state.counters["patches"]    = problem.patches;
        state.counters["age_groups"] = G;
        state.counters["steps"]      = integration_steps;
    }
    catch (const std::exception& exception) {
        state.SkipWithError(exception.what());
    }
}

void benchmark_cuda(benchmark::State& state)
{
    switch (static_cast<int>(state.range(1))) {
    case 1:
        benchmark_cuda_impl<1>(state);
        break;
    case 3:
        benchmark_cuda_impl<3>(state);
        break;
    case 6:
        benchmark_cuda_impl<6>(state);
        break;
    case 8:
        benchmark_cuda_impl<8>(state);
        break;
    default:
        state.SkipWithError("Unsupported age-group count.");
    }
}

void apply_shapes(benchmark::internal::Benchmark* benchmark)
{
    for (int patches : patch_counts) {
        for (int groups : age_group_counts) {
            benchmark->Args({patches, groups});
        }
    }
}

template <int G>
bool validate(std::string& error)
{
    constexpr int validation_patches = 263;
    ImplicitProblem cpu(validation_patches, G);
    ImplicitProblem gpu(validation_patches, G);
    advance_cpu<G>(cpu);
    GpuRunner<G> runner(gpu);
    runner.reset(gpu.initial_state);
    runner.run();
    const auto gpu_state = runner.download();

    double maximum_absolute_error = 0.0;
    double maximum_relative_error = 0.0;
    bool values_match             = true;
    for (size_t index = 0; index < cpu.state.size(); ++index) {
        const double absolute_error = std::abs(cpu.state[index] - gpu_state[index]);
        const double relative_error =
            absolute_error / std::max(validation_absolute_denominator, std::abs(cpu.state[index]));
        maximum_absolute_error = std::max(maximum_absolute_error, absolute_error);
        maximum_relative_error = std::max(maximum_relative_error, relative_error);
        if (!std::isfinite(gpu_state[index]) ||
            absolute_error > validation_relative_tolerance * (1.0 + std::abs(cpu.state[index]))) {
            values_match = false;
        }
    }
    for (int patch = 0; patch < validation_patches; ++patch) {
        for (int group = 0; group < G; ++group) {
            const size_t susceptible   = static_cast<size_t>(3 * group) * validation_patches + patch;
            const size_t exposed       = susceptible + validation_patches;
            const size_t infected      = exposed + validation_patches;
            const double recovered_cpu = cpu.population[cpu.group_patch_index(group, patch)] - cpu.state[susceptible] -
                                         cpu.state[exposed] - cpu.state[infected];
            const double recovered_gpu = gpu.population[gpu.group_patch_index(group, patch)] - gpu_state[susceptible] -
                                         gpu_state[exposed] - gpu_state[infected];
            const double absolute_error = std::abs(recovered_cpu - recovered_gpu);
            maximum_absolute_error      = std::max(maximum_absolute_error, absolute_error);
            maximum_relative_error =
                std::max(maximum_relative_error,
                         absolute_error / std::max(validation_absolute_denominator, std::abs(recovered_cpu)));
            if (!std::isfinite(recovered_gpu) ||
                absolute_error > validation_relative_tolerance * (1.0 + std::abs(recovered_cpu))) {
                values_match = false;
            }
        }
    }
    std::cout << "validation N_G=" << G << ": max_abs=" << maximum_absolute_error
              << ", max_rel=" << maximum_relative_error << '\n';
    if (!values_match) {
        error = "CPU and GPU results differ for N_G=" + std::to_string(G) +
                " (max_abs=" + std::to_string(maximum_absolute_error) +
                ", max_rel=" + std::to_string(maximum_relative_error) + ").";
    }
    return values_match;
}

bool validate_all(std::string& error)
{
    return validate<1>(error) && validate<3>(error) && validate<6>(error) && validate<8>(error);
}

} // namespace mio::benchmark_mio

BENCHMARK(mio::benchmark_mio::benchmark_serial)
    ->Apply(mio::benchmark_mio::apply_shapes)
    ->ArgNames({"patches", "age_groups"})
    ->Name("implicit/serial")
    ->UseRealTime();

BENCHMARK(mio::benchmark_mio::benchmark_cuda)
    ->Apply(mio::benchmark_mio::apply_shapes)
    ->ArgNames({"patches", "age_groups"})
    ->Name("implicit/cuda")
    ->UseRealTime();

int main(int argc, char** argv)
{
    bool needs_device = true;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (argument.find("--benchmark_list_tests") == 0 || argument == "--help" || argument == "-h") {
            needs_device = false;
        }
    }
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }
    if (needs_device) {
        try {
            cudaDeviceProp properties{};
            mio::benchmark_mio::check_cuda(cudaGetDeviceProperties(&properties, 0), "cudaGetDeviceProperties");
            std::cout << "GPU: " << properties.name << '\n';
            std::string error;
            if (!mio::benchmark_mio::validate_all(error)) {
                std::cerr << error << '\n';
                ::benchmark::Shutdown();
                return 1;
            }
        }
        catch (const std::exception& exception) {
            std::cerr << exception.what() << '\n';
            ::benchmark::Shutdown();
            return 1;
        }
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
    return 0;
}
