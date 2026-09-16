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
#include "ode_seir_runtime_explicit.h"
#include "ode_seir_roofline_cuda.h"

#include <cuda_runtime_api.h>

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

// stage_lambda is scratch space for patches * groups * 4 * temporal_block_steps doubles.
extern "C" cudaError_t launch_seir_stage_aligned_rk4_steps(double* totals, double* travelers, double* stage_lambda,
                                                           const double* contact_beta, const double* rate_exposed,
                                                           const double* rate_infected, int patches,
                                                           int travelers_per_patch, int groups, int steps, double dt,
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

template <class T, class Allocator>
cudaError_t upload(DeviceBuffer<T>& target, const std::vector<T, Allocator>& source, cudaStream_t stream = nullptr)
{
    return cudaMemcpyAsync(target.get(), source.data(), source.size() * sizeof(T), cudaMemcpyHostToDevice, stream);
}

template <class Expected, class Actual>
bool vectors_close(const Expected& expected, const Actual& actual)
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

bool validate_gpu_case(int groups, int patches, int travelers_per_patch, int steps, std::string& error)
{
    StageAlignedProblem reference(patches, travelers_per_patch, groups);
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
    DeviceBuffer<double> lambda(reference.stage_lambda.size() * temporal_block_steps);
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
    cudaStream_t stream        = nullptr;
    cudaGraph_t graph          = nullptr;
    cudaGraphExec_t executable = nullptr;
    const auto cleanup         = [&]() {
        if (executable) {
            cudaGraphExecDestroy(executable);
        }
        if (graph) {
            cudaGraphDestroy(graph);
        }
        if (stream) {
            cudaStreamDestroy(stream);
        }
    };
    auto status = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
    if (status == cudaSuccess) {
        status = cudaDeviceSynchronize();
    }
    if (status == cudaSuccess) {
        status = cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
    }
    if (status == cudaSuccess) {
        const auto launch_status = launch_seir_stage_aligned_rk4_steps(
            totals.get(), travelers.get(), lambda.get(), contact_beta.get(), rate_exposed.get(), rate_infected.get(),
            reference.patches, reference.travelers_per_patch, groups, steps, step_size, stream);
        const auto capture_status = cudaStreamEndCapture(stream, &graph);
        status                    = launch_status == cudaSuccess ? capture_status : launch_status;
    }
    if (status == cudaSuccess) {
        status = cudaGraphInstantiate(&executable, graph, 0);
    }
    if (status != cudaSuccess) {
        error = cuda_error(status, "CUDA validation graph setup failed");
        cleanup();
        return false;
    }
    auto expected = reference;
    advance_stage_aligned_reference(expected, step_size, steps);
    // Reuse the same graph after a reset: catches stale histories and reset errors.
    for (int repetition = 0; repetition < 2; ++repetition) {
        status = upload(totals, reference.totals, stream);
        if (status == cudaSuccess) {
            status = upload(travelers, reference.travelers, stream);
        }
        if (status == cudaSuccess) {
            status = cudaGraphLaunch(executable, stream);
        }
        if (status == cudaSuccess) {
            status = cudaStreamSynchronize(stream);
        }
        if (status != cudaSuccess) {
            error = cuda_error(status, "CUDA validation graph launch failed");
            cleanup();
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
            cleanup();
            return false;
        }

        if (!vectors_close(expected.totals, gpu_totals) || !vectors_close(expected.travelers, gpu_travelers)) {
            error = "Blocked CUDA graph and stepwise CPU RK4 differ for N_G=" + std::to_string(groups) +
                    ", patches=" + std::to_string(patches) + ", steps=" + std::to_string(steps);
            cleanup();
            return false;
        }
    }
    cleanup();
    return true;
}

bool validate_gpu(std::string& error)
{
    for (int groups : age_group_counts) {
        for (const auto& shape : {std::pair{3, 2}, std::pair{263, 262}}) {
            for (int steps : {1, 3, temporal_block_steps, temporal_block_steps + 1}) {
                if (!validate_gpu_case(groups, shape.first, shape.second, steps, error)) {
                    return false;
                }
            }
        }
        std::cout << "Stage-aligned CUDA graph vs stepwise CPU validation N_G=" << groups << ": passed\n";
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
        (problem.totals.size() + problem.travelers.size() + problem.stage_lambda.size() * temporal_block_steps +
         problem.contact_beta.size() + problem.rate_exposed.size() + problem.rate_infected.size()) *
        sizeof(double);
    if (status != cudaSuccess || required_memory > free_memory * 9 / 10) {
        state.SkipWithError("Insufficient device memory for this benchmark shape.");
        return;
    }

    DeviceBuffer<double> totals(problem.totals.size());
    DeviceBuffer<double> travelers(problem.travelers.size());
    DeviceBuffer<double> lambda(problem.stage_lambda.size() * temporal_block_steps);
    DeviceBuffer<double> contact_beta(problem.contact_beta.size());
    DeviceBuffer<double> rate_exposed(problem.rate_exposed.size());
    DeviceBuffer<double> rate_infected(problem.rate_infected.size());
    if (!allocate_successfully(state, {totals.status(), travelers.status(), lambda.status(), contact_beta.status(),
                                       rate_exposed.status(), rate_infected.status()})) {
        return;
    }

    cudaStream_t stream             = nullptr;
    cudaGraph_t graph               = nullptr;
    cudaGraphExec_t graph_exec      = nullptr;
    const auto destroy_cuda_objects = [&]() {
        if (graph_exec != nullptr) {
            cudaGraphExecDestroy(graph_exec);
        }
        if (graph != nullptr) {
            cudaGraphDestroy(graph);
        }
        if (stream != nullptr) {
            cudaStreamDestroy(stream);
        }
    };

    status = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
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
        destroy_cuda_objects();
        return;
    }

    const auto enqueue_steps = [&](int steps) {
        return launch_seir_stage_aligned_rk4_steps(
            totals.get(), travelers.get(), lambda.get(), contact_beta.get(), rate_exposed.get(), rate_infected.get(),
            problem.patches, problem.travelers_per_patch, problem.groups, steps, step_size, stream);
    };

    status = upload(totals, problem.totals, stream);
    if (status == cudaSuccess) {
        status = upload(travelers, problem.travelers, stream);
    }
    if (status == cudaSuccess) {
        status = enqueue_steps(1);
    }
    if (status == cudaSuccess) {
        status = cudaStreamSynchronize(stream);
    }
    if (status == cudaSuccess) {
        status = cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
    }
    if (status == cudaSuccess) {
        const auto capture_status = enqueue_steps(integration_steps);
        const auto end_status     = cudaStreamEndCapture(stream, &graph);
        status                    = capture_status == cudaSuccess ? end_status : capture_status;
    }
    if (status == cudaSuccess) {
        status = cudaGraphInstantiate(&graph_exec, graph, 0);
    }
    if (status == cudaSuccess) {
        status = cudaGraphUpload(graph_exec, stream);
    }
    if (status == cudaSuccess) {
        status = cudaStreamSynchronize(stream);
    }
    if (status != cudaSuccess) {
        const auto message = cuda_error(status, "CUDA graph setup failed");
        state.SkipWithError(message.c_str());
        destroy_cuda_objects();
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
            status = cudaGraphLaunch(graph_exec, stream);
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
    destroy_cuda_objects();

    state.counters["patches"]                = problem.patches;
    state.counters["edges"]                  = static_cast<double>(problem.edges());
    state.counters["age_groups"]             = problem.groups;
    state.counters["steps"]                  = integration_steps;
    state.counters["implementation_version"] = stage_aligned_implementation_version;
    state.counters["temporal_block_steps"]   = temporal_block_steps;
    state.SetItemsProcessed(state.iterations() * integration_steps * static_cast<int64_t>(problem.edges()));
}

extern "C" cudaError_t launch_seir_daily_home(double*, const double*, const double*, const double*, int, int, int,
                                              double, cudaStream_t);
extern "C" cudaError_t launch_seir_daily_departure(double*, double*, double*, const double*, int, int, cudaStream_t);
extern "C" cudaError_t launch_seir_daily_return(double*, const double*, double*, int, int, cudaStream_t);
namespace scenario = mio::runtime_scenario;

class RuntimeGpuRunner
{
public:
    RuntimeGpuRunner(const scenario::Inputs& inputs, scenario::Schedule schedule)
        : patches(inputs.p)
        , groups(inputs.g)
        , time(schedule)
        , totals(inputs.initial.size())
        , travelers(static_cast<size_t>(inputs.p) * (inputs.p - 1) * 4 * inputs.g)
        , history(inputs.initial.size() * temporal_block_steps)
        , snapshot(inputs.initial.size())
        , ht(inputs.ht.size())
        , beta(inputs.beta.size())
        , rate_e(inputs.rate_e.size())
        , rate_i(inputs.rate_i.size())
    {
        try {
            for (auto status : {totals.status(), travelers.status(), history.status(), snapshot.status(), ht.status(),
                                beta.status(), rate_e.status(), rate_i.status()})
                require(status, "allocate daily scenario");
            require(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking), "create daily stream");
            require(upload(ht, inputs.ht, stream), "upload H");
            require(upload(beta, inputs.beta, stream), "upload beta");
            require(upload(rate_e, inputs.rate_e, stream), "upload exposed rates");
            require(upload(rate_i, inputs.rate_i, stream), "upload infected rates");
            reset(inputs.initial);
            enqueue_day();
            require(cudaStreamSynchronize(stream), "warm up daily scenario");
            reset(inputs.initial);
            require(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal), "begin daily graph");
            try {
                enqueue_day();
            }
            catch (...) {
                cudaGraph_t failed = nullptr;
                cudaStreamEndCapture(stream, &failed);
                if (failed)
                    cudaGraphDestroy(failed);
                throw;
            }
            require(cudaStreamEndCapture(stream, &graph), "end daily graph");
            require(cudaGraphInstantiate(&executable, graph, 0), "instantiate daily graph");
            require(cudaGraphUpload(executable, stream), "upload daily graph");
            require(cudaStreamSynchronize(stream), "finish daily setup");
        }
        catch (...) {
            release();
            throw;
        }
    }
    ~RuntimeGpuRunner()
    {
        release();
    }
    RuntimeGpuRunner(const RuntimeGpuRunner&) = delete;
    RuntimeGpuRunner& operator=(const RuntimeGpuRunner&) = delete;
    void reset(const std::vector<double>& initial)
    {
        require(upload(totals, initial, stream), "reset daily residents");
        require(cudaStreamSynchronize(stream), "finish daily reset");
    }
    void run(bool profile = false)
    {
#ifdef MEMILIO_BENCHMARK_ROOFLINE
        if (profile) {
            mio::benchmark_roofline::run_profiled_cuda_days(executable, stream, time.days);
            return;
        }
#else
        (void)profile;
#endif
        for (int day = 0; day < time.days; ++day)
            require(cudaGraphLaunch(executable, stream), "launch daily graph");
        require(cudaStreamSynchronize(stream), "finish daily simulation");
    }
    std::vector<double> download() const
    {
        std::vector<double> result(static_cast<size_t>(4 * groups) * patches);
        require(cudaMemcpy(result.data(), totals.get(), result.size() * sizeof(double), cudaMemcpyDeviceToHost),
                "download daily residents");
        return result;
    }
    const double* device_state() const
    {
        return totals.get();
    }

private:
    static void require(cudaError_t status, const char* operation)
    {
        if (status != cudaSuccess)
            throw std::runtime_error(cuda_error(status, operation));
    }
    void enqueue_day()
    {
        require(launch_seir_daily_home(totals.get(), beta.get(), rate_e.get(), rate_i.get(), patches, groups,
                                       time.half_day_steps, time.dt(), stream),
                "integrate home phase");
        require(launch_seir_daily_departure(totals.get(), travelers.get(), snapshot.get(), ht.get(), patches, groups,
                                            stream),
                "depart");
        require(launch_seir_stage_aligned_rk4_steps(totals.get(), travelers.get(), history.get(), beta.get(),
                                                    rate_e.get(), rate_i.get(), patches, patches - 1, groups,
                                                    time.half_day_steps, time.dt(), stream),
                "integrate away phase");
        require(launch_seir_daily_return(totals.get(), travelers.get(), snapshot.get(), patches, groups, stream),
                "return home");
    }
    void release() noexcept
    {
        if (executable)
            cudaGraphExecDestroy(executable);
        if (graph)
            cudaGraphDestroy(graph);
        if (stream)
            cudaStreamDestroy(stream);
    }
    int patches, groups;
    scenario::Schedule time;
    DeviceBuffer<double> totals, travelers, history, snapshot, ht, beta, rate_e, rate_i;
    cudaStream_t stream        = nullptr;
    cudaGraph_t graph          = nullptr;
    cudaGraphExec_t executable = nullptr;
};

void runtime_cuda(benchmark::State& state)
{
    try {
        const bool profile = mio::benchmark_roofline::enabled();
        if (profile)
            (void)mio::benchmark_roofline::selected_day();
        scenario::Inputs inputs(static_cast<int>(state.range(0)), static_cast<int>(state.range(1)));
        RuntimeGpuRunner runner(inputs, scenario::schedule());
        bool profile_window_completed = false;
        for (auto _ : state) {
            state.PauseTiming();
            runner.reset(inputs.initial);
            state.ResumeTiming();
            runner.run(profile);
            profile_window_completed = profile;
            benchmark::DoNotOptimize(runner.device_state());
        }
        scenario::check_population(inputs, runner.download());
        scenario::counters(state, inputs, true, 0);
        if (profile)
            state.counters["roofline_profile_window_completed"] = profile_window_completed ? 1 : 0;
    }
    catch (const std::exception& error) {
        state.SkipWithError(error.what());
    }
}

template <int G>
void validate_runtime_gpu()
{
    for (int patches : {3, 263}) {
        for (int phase_steps : {1, 32, 65}) {
            for (bool no_mobility : {false, true}) {
                if (patches == 263 && no_mobility)
                    continue;
                const scenario::Inputs inputs(patches, G, no_mobility);
                const scenario::Schedule time{2, phase_steps};
                scenario::ExplicitProblem expected(inputs);
                scenario::advance_explicit<G>(expected, time);
                if (patches == 3)
                    scenario::compare(inputs, scenario::Reference(inputs, true).run(time), expected.core.totals);
                RuntimeGpuRunner runner(inputs, time);
                for (int replay = 0; replay < 2; ++replay) {
                    runner.reset(inputs.initial);
                    runner.run();
                    scenario::compare(inputs, expected.core.totals, runner.download());
                }
            }
        }
    }
    std::cout << "Runtime explicit daily CUDA graph N_G=" << G << ": passed\n";
}

void apply_cuda_shapes(benchmark::internal::Benchmark* benchmark)
{
    for (const auto& [patches, travelers] : problem_shapes) {
        for (int groups : age_group_counts) {
            benchmark->Args({patches, travelers, groups});
        }
    }
    benchmark->Args({stage_aligned_strong_scaling_patches, stage_aligned_strong_scaling_patches - 1, 6});
}

} // namespace mio::benchmark_mio

BENCHMARK(mio::benchmark_mio::benchmark_stage_aligned_cuda)
    ->Apply(mio::benchmark_mio::apply_cuda_shapes)
    ->ArgNames({"patches", "travelers_per_patch", "age_groups"})
    ->Name("stage_aligned/cuda")
    ->UseRealTime();

int main(int argc, char** argv)
{
    const bool runtime = mio::runtime_scenario::enabled();
    if (runtime)
        mio::runtime_scenario::register_shapes("runtime/explicit/cuda", mio::benchmark_mio::runtime_cuda);
    bool needs_device = true;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (mio::runtime_scenario::informational_argument(argument)) {
            needs_device = false;
        }
    }
    ::benchmark::Initialize(&argc, argv);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }
    std::string error;
    try {
        if (runtime && needs_device) {
            if (mio::benchmark_roofline::enabled())
                (void)mio::benchmark_roofline::selected_day();
            mio::runtime_scenario::validate_accuracy();
            mio::benchmark_mio::validate_runtime_gpu<1>();
            mio::benchmark_mio::validate_runtime_gpu<3>();
            mio::benchmark_mio::validate_runtime_gpu<6>();
            mio::benchmark_mio::validate_runtime_gpu<8>();
        }
    }
    catch (const std::exception& exception) {
        std::cerr << "Runtime CUDA validation failed: " << exception.what() << '\n';
        ::benchmark::Shutdown();
        return 1;
    }
    if (!runtime && needs_device && !mio::benchmark_mio::validate_gpu(error)) {
        std::cerr << error << '\n';
        ::benchmark::Shutdown();
        return 1;
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();
    return 0;
}
