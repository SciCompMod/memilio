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
#ifndef MIO_ODE_SEIR_ROOFLINE_CUDA_H
#define MIO_ODE_SEIR_ROOFLINE_CUDA_H

#include "ode_seir_roofline.h"

#ifdef MEMILIO_BENCHMARK_ROOFLINE
#include <cuda_profiler_api.h>
#include <cuda_runtime_api.h>
#include <stdexcept>
#include <string>

namespace mio::benchmark_roofline
{
inline void check_cuda_profile(cudaError_t status, const char* operation)
{
    if (status != cudaSuccess)
        throw std::runtime_error(std::string("Roofline ") + operation + ": " + cudaGetErrorString(status));
}

inline void run_profiled_cuda_days(cudaGraphExec_t graph, cudaStream_t stream, int days)
{
    const int selected = selected_day();
    if (selected > days)
        throw std::invalid_argument("ROOFLINE_DAY exceeds the number of daily CUDA graph replays.");
    for (int day = 1; day <= days; ++day) {
        if (day != selected) {
            check_cuda_profile(cudaGraphLaunch(graph, stream), "launch unprofiled daily graph");
            continue;
        }
        // The capture and warm-up are already complete. Drain preceding days
        // before opening the profiler window; keep the captured graph intact.
        check_cuda_profile(cudaStreamSynchronize(stream), "finish preceding daily graphs");
        check_cuda_profile(cudaProfilerStart(), "start selected-day profiling");
        try {
            check_cuda_profile(cudaGraphLaunch(graph, stream), "launch selected daily graph");
            check_cuda_profile(cudaStreamSynchronize(stream), "finish selected daily graph");
        }
        catch (...) {
            // Always attempt to close an opened window, including asynchronous
            // CUDA failures. Preserve the original launch/synchronization error.
            (void)cudaProfilerStop();
            throw;
        }
        check_cuda_profile(cudaProfilerStop(), "stop selected-day profiling");
    }
    check_cuda_profile(cudaStreamSynchronize(stream), "finish remaining daily graphs");
}
} // namespace mio::benchmark_roofline
#endif

#endif
