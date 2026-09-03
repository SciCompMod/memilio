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
#include "ode_seir_benchmark_stage_aligned.h"

#include <cuda_runtime_api.h>

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

extern "C" cudaError_t launch_seir_stage_aligned_rk4_step(double* totals, double* travelers, double* stage_lambda,
                                                          const double* contact_beta, const double* rate_exposed,
                                                          const double* rate_infected, int patches,
                                                          int travelers_per_patch, int groups, double dt,
                                                          cudaStream_t stream);

namespace mio::benchmark_mio
{

template <class T>
class DeviceBuffer
{
public:
    explicit DeviceBuffer(size_t count)
        : m_status(cudaMalloc(reinterpret_cast<void**>(&m_data), count * sizeof(T)))
    {
    }

    ~DeviceBuffer()
    {
        cudaFree(m_data);
    }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    T* get()
    {
        return m_data;
    }

    const T* get() const
    {
        return m_data;
    }

    cudaError_t status() const
    {
        return m_status;
    }

private:
    T* m_data = nullptr;
    cudaError_t m_status;
};

std::string cuda_error(cudaError_t status, const char* operation)
{
    return std::string(operation) + ": " + cudaGetErrorString(status);
}

bool allocate_successfully(benchmark::State& state, std::initializer_list<cudaError_t> statuses)
{
    for (auto status : statuses) {
        if (status != cudaSuccess) {
            const auto message = cuda_error(status, "CUDA allocation failed");
            state.SkipWithError(message.c_str());
            return false;
        }
    }
    return true;
}

template <class T>
cudaError_t upload(DeviceBuffer<T>& target, const std::vector<T>& source, cudaStream_t stream = nullptr)
{
    return cudaMemcpyAsync(target.get(), source.data(), source.size() * sizeof(T), cudaMemcpyHostToDevice, stream);
}

bool vectors_close(const std::vector<double>& expected, const std::vector<double>& actual)
{
    if (expected.size() != actual.size()) {
        return false;
    }
    for (size_t index = 0; index < expected.size(); ++index) {
        if (!std::isfinite(expected[index]) || !std::isfinite(actual[index]) ||
            std::abs(expected[index] - actual[index]) > 1e-10 * (1.0 + std::abs(expected[index]))) {
            return false;
        }
    }
    return true;
}

bool validate_gpu(std::string& error)
{
    for (int groups : age_group_counts) {
        StageAlignedProblem reference(3, 2, groups);
        for (int compartment = 0; compartment < reference.compartments(); ++compartment) {
            for (int patch = 0; patch < reference.patches; ++patch) {
                reference.totals[static_cast<size_t>(compartment) * reference.patches + patch] *=
                    1.0 + 0.01 * ((3 * compartment + patch) % 7);
                for (int traveler = 0; traveler < reference.travelers_per_patch; ++traveler) {
                    const size_t index = (static_cast<size_t>(patch) * reference.compartments() + compartment) *
                                             reference.travelers_per_patch +
                                         traveler;
                    reference.travelers[index] *= 0.9 + 0.01 * ((compartment + 2 * patch + traveler) % 9);
                }
            }
        }
        for (int group = 0; group < groups; ++group) {
            reference.rate_exposed[group] *= 1.0 + 0.02 * group;
            reference.rate_infected[group] *= 1.0 + 0.03 * group;
        }
        DeviceBuffer<double> totals(reference.totals.size());
        DeviceBuffer<double> travelers(reference.travelers.size());
        DeviceBuffer<double> lambda(reference.stage_lambda.size());
        DeviceBuffer<double> contact_beta(reference.contact_beta.size());
        DeviceBuffer<double> rate_exposed(reference.rate_exposed.size());
        DeviceBuffer<double> rate_infected(reference.rate_infected.size());
        for (auto status : {totals.status(), travelers.status(), lambda.status(), contact_beta.status(),
                            rate_exposed.status(), rate_infected.status()}) {
            if (status != cudaSuccess) {
                error = cuda_error(status, "CUDA validation allocation failed");
                return false;
            }
        }

        for (auto status : {upload(totals, reference.totals), upload(travelers, reference.travelers),
                            upload(contact_beta, reference.contact_beta), upload(rate_exposed, reference.rate_exposed),
                            upload(rate_infected, reference.rate_infected)}) {
            if (status != cudaSuccess) {
                error = cuda_error(status, "CUDA validation upload failed");
                return false;
            }
        }
        auto status = launch_seir_stage_aligned_rk4_step(
            totals.get(), travelers.get(), lambda.get(), contact_beta.get(), rate_exposed.get(), rate_infected.get(),
            reference.patches, reference.travelers_per_patch, groups, step_size, nullptr);
        if (status == cudaSuccess) {
            status = cudaDeviceSynchronize();
        }
        if (status != cudaSuccess) {
            error = cuda_error(status, "CUDA validation kernel failed");
            return false;
        }

        std::vector<double> gpu_totals(reference.totals.size());
        std::vector<double> gpu_travelers(reference.travelers.size());
        status =
            cudaMemcpy(gpu_totals.data(), totals.get(), gpu_totals.size() * sizeof(double), cudaMemcpyDeviceToHost);
        if (status == cudaSuccess) {
            status = cudaMemcpy(gpu_travelers.data(), travelers.get(), gpu_travelers.size() * sizeof(double),
                                cudaMemcpyDeviceToHost);
        }
        if (status != cudaSuccess) {
            error = cuda_error(status, "CUDA validation download failed");
            return false;
        }

        advance_stage_aligned(reference);
        if (!vectors_close(reference.totals, gpu_totals) || !vectors_close(reference.travelers, gpu_travelers)) {
            error = "CUDA and CPU stage-aligned RK4 results differ.";
            return false;
        }
    }
    return true;
}

void benchmark_stage_aligned_cuda(benchmark::State& state)
{
    StageAlignedProblem problem(static_cast<int>(state.range(0)), static_cast<int>(state.range(1)),
                                static_cast<int>(state.range(2)));
    size_t free_memory  = 0;
    size_t total_memory = 0;
    auto status         = cudaMemGetInfo(&free_memory, &total_memory);
    const size_t required_memory =
        (problem.totals.size() + problem.travelers.size() + problem.stage_lambda.size() + problem.contact_beta.size() +
         problem.rate_exposed.size() + problem.rate_infected.size()) *
        sizeof(double);
    if (status != cudaSuccess || required_memory > free_memory * 9 / 10) {
        state.SkipWithError("Insufficient device memory for this benchmark shape.");
        return;
    }

    DeviceBuffer<double> totals(problem.totals.size());
    DeviceBuffer<double> travelers(problem.travelers.size());
    DeviceBuffer<double> lambda(problem.stage_lambda.size());
    DeviceBuffer<double> contact_beta(problem.contact_beta.size());
    DeviceBuffer<double> rate_exposed(problem.rate_exposed.size());
    DeviceBuffer<double> rate_infected(problem.rate_infected.size());
    if (!allocate_successfully(state, {totals.status(), travelers.status(), lambda.status(), contact_beta.status(),
                                       rate_exposed.status(), rate_infected.status()})) {
        return;
    }

    cudaStream_t stream = nullptr;
    status              = cudaStreamCreate(&stream);
    if (status == cudaSuccess) {
        status = upload(contact_beta, problem.contact_beta, stream);
    }
    if (status == cudaSuccess) {
        status = upload(rate_exposed, problem.rate_exposed, stream);
    }
    if (status == cudaSuccess) {
        status = upload(rate_infected, problem.rate_infected, stream);
    }
    if (status == cudaSuccess) {
        status = cudaStreamSynchronize(stream);
    }
    if (status != cudaSuccess) {
        const auto message = cuda_error(status, "CUDA setup failed");
        state.SkipWithError(message.c_str());
        if (stream != nullptr) {
            cudaStreamDestroy(stream);
        }
        return;
    }

    for (auto _ : state) {
        state.PauseTiming();
        problem.reset_state();
        status = upload(totals, problem.totals, stream);
        if (status == cudaSuccess) {
            status = upload(travelers, problem.travelers, stream);
        }
        if (status == cudaSuccess) {
            status = cudaStreamSynchronize(stream);
        }
        state.ResumeTiming();

        if (status == cudaSuccess) {
            for (int step = 0; step < integration_steps && status == cudaSuccess; ++step) {
                status =
                    launch_seir_stage_aligned_rk4_step(totals.get(), travelers.get(), lambda.get(), contact_beta.get(),
                                                       rate_exposed.get(), rate_infected.get(), problem.patches,
                                                       problem.travelers_per_patch, problem.groups, step_size, stream);
            }
        }
        if (status == cudaSuccess) {
            status = cudaStreamSynchronize(stream);
        }
        if (status != cudaSuccess) {
            const auto message = cuda_error(status, "CUDA benchmark step failed");
            state.SkipWithError(message.c_str());
            break;
        }
        benchmark::DoNotOptimize(travelers.get());
    }
    if (stream != nullptr) {
        cudaStreamDestroy(stream);
    }

    state.counters["patches"]    = problem.patches;
    state.counters["edges"]      = static_cast<double>(problem.edges());
    state.counters["age_groups"] = problem.groups;
    state.counters["steps"]      = integration_steps;
    state.SetItemsProcessed(state.iterations() * integration_steps * static_cast<int64_t>(problem.edges()));
}

void apply_cuda_shapes(benchmark::internal::Benchmark* benchmark)
{
    for (const auto& [patches, travelers] : problem_shapes) {
        for (int groups : age_group_counts) {
            benchmark->Args({patches, travelers, groups});
        }
    }
}

} // namespace mio::benchmark_mio

BENCHMARK(mio::benchmark_mio::benchmark_stage_aligned_cuda)
    ->Apply(mio::benchmark_mio::apply_cuda_shapes)
    ->ArgNames({"patches", "travelers_per_patch", "age_groups"})
    ->Name("stage_aligned/cuda")
    ->UseRealTime();

int main(int argc, char** argv)
{
    bool needs_device = true;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (argument.starts_with("--benchmark_list_tests") || argument == "--help" || argument == "-h") {
            needs_device = false;
        }
    }
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }
    std::string error;
    if (needs_device && !mio::benchmark_mio::validate_gpu(error)) {
        std::cerr << error << '\n';
        ::benchmark::Shutdown();
        return 1;
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
    return 0;
}
