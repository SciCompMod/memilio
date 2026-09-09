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

#include <cuda_runtime.h>
#include "ode_seir_benchmark_stage_aligned.h"

#include <cstddef>

namespace
{

constexpr int block_size = 256;
constexpr int time_tile  = mio::benchmark_mio::temporal_block_steps;

template <int G, bool StoreHistory = true>
__global__ void
seir_totals_allpatches_ct_kernel(double* __restrict__ totals, double* __restrict__ stage_lambda,
                                 const double* __restrict__ contact_beta, const double* __restrict__ rate_exposed,
                                 const double* __restrict__ rate_infected, int patches, int steps, double dt)
{
    __shared__ double parameters[G * G + 2 * G];
    for (int index = static_cast<int>(threadIdx.x); index < G * G + 2 * G; index += static_cast<int>(blockDim.x)) {
        if (index < G * G) {
            parameters[index] = contact_beta[index];
        }
        else if (index < G * G + G) {
            parameters[index] = rate_exposed[index - G * G];
        }
        else {
            parameters[index] = rate_infected[index - G * G - G];
        }
    }
    __syncthreads();

    const int patch = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (patch >= patches) {
        return;
    }

    // Aggregate trajectories do not depend on the traveler states. Produce the
    // complete block of RK4 stage rates before launching the traveler kernel.
    for (int step = 0; step < steps; ++step) {
        double s0[G], e0[G], i0[G];
        double s[G], e[G], i[G];
        double result_s[G], result_e[G], result_i[G], inverse_population[G], lambda[G];
#pragma unroll
        for (int group = 0; group < G; ++group) {
            s0[group] = s[group] = result_s[group] = totals[(4 * group) * patches + patch];
            e0[group] = e[group] = result_e[group] = totals[(4 * group + 1) * patches + patch];
            i0[group] = i[group] = result_i[group] = totals[(4 * group + 2) * patches + patch];
            const double population   = s0[group] + e0[group] + i0[group] + totals[(4 * group + 3) * patches + patch];
            inverse_population[group] = population > 1e-12 ? 1.0 / population : 0.0;
        }

        const double sixth_dt = dt / 6.0;
        const double third_dt = dt / 3.0;
        const double half_dt  = 0.5 * dt;
        // Keep the stage loop rolled to reduce register pressure and spilling
        // from fully unrolling all stages and groups simultaneously.
#pragma unroll 1
        for (int stage = 0; stage < 4; ++stage) {
#pragma unroll
            for (int target = 0; target < G; ++target) {
                lambda[target] = 0.0;
#pragma unroll
                for (int source = 0; source < G; ++source) {
                    lambda[target] =
                        fma(parameters[target * G + source], i[source] * inverse_population[source], lambda[target]);
                }
                if constexpr (StoreHistory) {
                    stage_lambda[((static_cast<size_t>(patch) * G + target) * time_tile + step) * 4 + stage] =
                        lambda[target];
                }
            }

            const double weight     = (stage == 0 || stage == 3) ? sixth_dt : third_dt;
            const double stage_step = stage < 2 ? half_dt : dt;
#pragma unroll
            for (int group = 0; group < G; ++group) {
                const double flow_se = lambda[group] * s[group];
                const double flow_ei = parameters[G * G + group] * e[group];
                const double flow_ir = parameters[G * G + G + group] * i[group];
                const double ds      = -flow_se;
                const double de      = flow_se - flow_ei;
                const double di      = flow_ei - flow_ir;
                result_s[group]      = fma(weight, ds, result_s[group]);
                result_e[group]      = fma(weight, de, result_e[group]);
                result_i[group]      = fma(weight, di, result_i[group]);
                if (stage < 3) {
                    s[group] = fma(stage_step, ds, s0[group]);
                    e[group] = fma(stage_step, de, e0[group]);
                    i[group] = fma(stage_step, di, i0[group]);
                }
            }
        }

#pragma unroll
        for (int group = 0; group < G; ++group) {
            const size_t susceptible = static_cast<size_t>(4 * group) * patches + patch;
            const size_t exposed     = static_cast<size_t>(4 * group + 1) * patches + patch;
            const size_t infected    = static_cast<size_t>(4 * group + 2) * patches + patch;
            const size_t recovered   = static_cast<size_t>(4 * group + 3) * patches + patch;
            const double population  = totals[susceptible] + totals[exposed] + totals[infected] + totals[recovered];
            totals[susceptible]      = result_s[group];
            totals[exposed]          = result_e[group];
            totals[infected]         = result_i[group];
            totals[recovered]        = population - result_s[group] - result_e[group] - result_i[group];
        }
    }
}

__device__ __forceinline__ void integrate_traveler_group(double& susceptible, double& exposed, double& infected,
                                                         double& recovered, const double* lambda, double rate_exposed,
                                                         double rate_infected, double dt)
{
    const double s0         = susceptible;
    const double e0         = exposed;
    const double i0         = infected;
    const double population = s0 + e0 + i0 + recovered;
    const double half_dt    = 0.5 * dt;
    const double sixth_dt   = dt / 6.0;
    const double third_dt   = dt / 3.0;

    double flow_se  = lambda[0] * s0;
    double flow_ei  = rate_exposed * e0;
    double flow_ir  = rate_infected * i0;
    double result_s = fma(-sixth_dt, flow_se, s0);
    double result_e = fma(sixth_dt, flow_se - flow_ei, e0);
    double result_i = fma(sixth_dt, flow_ei - flow_ir, i0);
    double s        = fma(-half_dt, flow_se, s0);
    double e        = fma(half_dt, flow_se - flow_ei, e0);
    double i        = fma(half_dt, flow_ei - flow_ir, i0);

    flow_se  = lambda[1] * s;
    flow_ei  = rate_exposed * e;
    flow_ir  = rate_infected * i;
    result_s = fma(-third_dt, flow_se, result_s);
    result_e = fma(third_dt, flow_se - flow_ei, result_e);
    result_i = fma(third_dt, flow_ei - flow_ir, result_i);
    s        = fma(-half_dt, flow_se, s0);
    e        = fma(half_dt, flow_se - flow_ei, e0);
    i        = fma(half_dt, flow_ei - flow_ir, i0);

    flow_se  = lambda[2] * s;
    flow_ei  = rate_exposed * e;
    flow_ir  = rate_infected * i;
    result_s = fma(-third_dt, flow_se, result_s);
    result_e = fma(third_dt, flow_se - flow_ei, result_e);
    result_i = fma(third_dt, flow_ei - flow_ir, result_i);
    s        = fma(-dt, flow_se, s0);
    e        = fma(dt, flow_se - flow_ei, e0);
    i        = fma(dt, flow_ei - flow_ir, i0);

    flow_se     = lambda[3] * s;
    flow_ei     = rate_exposed * e;
    flow_ir     = rate_infected * i;
    susceptible = fma(-sixth_dt, flow_se, result_s);
    exposed     = fma(sixth_dt, flow_se - flow_ei, result_e);
    infected    = fma(sixth_dt, flow_ei - flow_ir, result_i);
    recovered   = population - susceptible - exposed - infected;
}

template <int G>
__global__ void
seir_traveler_allpatches_ct_kernel(double* __restrict__ travelers, const double* __restrict__ stage_lambda,
                                   const double* __restrict__ rate_exposed, const double* __restrict__ rate_infected,
                                   int travelers_per_patch, int steps, double dt)
{
    // One block owns a patch/group/traveler tile. Shared rates are contiguous;
    // the four traveler states remain in registers for the whole time block.
    __shared__ double parameters[4 * time_tile + 2];
    const int patch             = static_cast<int>(blockIdx.y);
    const int group             = static_cast<int>(blockIdx.z);
    const size_t history_offset = (static_cast<size_t>(patch) * G + group) * time_tile * 4;
    for (int index = static_cast<int>(threadIdx.x); index < 4 * steps; index += static_cast<int>(blockDim.x)) {
        parameters[index] = stage_lambda[history_offset + index];
    }
    if (threadIdx.x == 0) {
        parameters[4 * time_tile]     = rate_exposed[group];
        parameters[4 * time_tile + 1] = rate_infected[group];
    }
    __syncthreads();

    const int traveler = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (traveler >= travelers_per_patch) {
        return;
    }

    const size_t patch_offset      = static_cast<size_t>(patch) * 4 * G * travelers_per_patch;
    const size_t susceptible_index = patch_offset + static_cast<size_t>(4 * group) * travelers_per_patch + traveler;
    const size_t exposed_index     = susceptible_index + travelers_per_patch;
    const size_t infected_index    = exposed_index + travelers_per_patch;
    const size_t recovered_index   = infected_index + travelers_per_patch;
    double susceptible             = travelers[susceptible_index];
    double exposed                 = travelers[exposed_index];
    double infected                = travelers[infected_index];
    double recovered               = travelers[recovered_index];
    for (int step = 0; step < steps; ++step) {
        integrate_traveler_group(susceptible, exposed, infected, recovered, parameters + 4 * step,
                                 parameters[4 * time_tile], parameters[4 * time_tile + 1], dt);
    }
    travelers[susceptible_index] = susceptible;
    travelers[exposed_index]     = exposed;
    travelers[infected_index]    = infected;
    travelers[recovered_index]   = recovered;
}

template <int G>
cudaError_t launch_step(double* totals, double* travelers, double* stage_lambda, const double* contact_beta,
                        const double* rate_exposed, const double* rate_infected, int patches, int travelers_per_patch,
                        int steps, double dt, cudaStream_t stream)
{
    for (int completed = 0; completed < steps;) {
        const int count    = std::min(time_tile, steps - completed);
        int totals_threads = 32;
        // Aggregate work is patch-parallel only. Small blocks distribute the
        // relatively few patch trajectories over more SMs than 256-thread blocks.
        while (totals_threads < patches && totals_threads < 64) {
            totals_threads *= 2;
        }
        const int totals_blocks = (patches + totals_threads - 1) / totals_threads;
        seir_totals_allpatches_ct_kernel<G><<<totals_blocks, totals_threads, 0, stream>>>(
            totals, stage_lambda, contact_beta, rate_exposed, rate_infected, patches, count, dt);
        auto status = cudaGetLastError();
        if (status != cudaSuccess) {
            return status;
        }

        int traveler_threads = 32;
        while (traveler_threads < travelers_per_patch && traveler_threads < block_size) {
            traveler_threads *= 2;
        }
        const dim3 traveler_grid(
            static_cast<unsigned int>((travelers_per_patch + traveler_threads - 1) / traveler_threads),
            static_cast<unsigned int>(patches), G);
        seir_traveler_allpatches_ct_kernel<G><<<traveler_grid, traveler_threads, 0, stream>>>(
            travelers, stage_lambda, rate_exposed, rate_infected, travelers_per_patch, count, dt);
        status = cudaGetLastError();
        if (status != cudaSuccess) {
            return status;
        }
        completed += count;
    }
    return cudaSuccess;
}

// The following kernels are used only by the daily runtime scenario. H is
// destination-major, while resident/totals arrays are compartment-major.
__global__ void daily_departure_kernel(const double* residents, const double* ht, double* totals, double* travelers,
                                       int patches, int compartments)
{
    const int dest = static_cast<int>(blockIdx.x), c = static_cast<int>(blockIdx.y);
    const int count = patches - 1;
    auto* row       = travelers + (static_cast<size_t>(dest) * compartments + c) * count;
    double sum      = 0.0;
    for (int t = static_cast<int>(threadIdx.x); t < count; t += static_cast<int>(blockDim.x)) {
        const int origin = t + (t >= dest);
        const double value =
            ht[static_cast<size_t>(dest) * patches + origin] * residents[static_cast<size_t>(c) * patches + origin];
        row[t] = value;
        sum += value;
    }
    __shared__ double partial[256];
    partial[threadIdx.x] = sum;
    __syncthreads();
    for (int stride = static_cast<int>(blockDim.x) / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < static_cast<unsigned>(stride))
            partial[threadIdx.x] += partial[threadIdx.x + stride];
        __syncthreads();
    }
    if (threadIdx.x == 0)
        totals[static_cast<size_t>(c) * patches + dest] =
            partial[0] +
            ht[static_cast<size_t>(dest) * patches + dest] * residents[static_cast<size_t>(c) * patches + dest];
}

__global__ void daily_stayers_kernel(const double* totals, const double* travelers, double* stayers, int patches,
                                     int compartments)
{
    const int dest = static_cast<int>(blockIdx.x), c = static_cast<int>(blockIdx.y);
    const auto* row = travelers + (static_cast<size_t>(dest) * compartments + c) * (patches - 1);
    double sum      = 0.0;
    for (int t = static_cast<int>(threadIdx.x); t < patches - 1; t += static_cast<int>(blockDim.x))
        sum += row[t];
    __shared__ double partial[256];
    partial[threadIdx.x] = sum;
    __syncthreads();
    for (int stride = static_cast<int>(blockDim.x) / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < static_cast<unsigned>(stride))
            partial[threadIdx.x] += partial[threadIdx.x + stride];
        __syncthreads();
    }
    if (threadIdx.x == 0)
        stayers[static_cast<size_t>(c) * patches + dest] = totals[static_cast<size_t>(c) * patches + dest] - partial[0];
}

__global__ void daily_return_kernel(const double* stayers, const double* travelers, double* residents, int patches,
                                    int compartments)
{
    const int origin = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x), c = static_cast<int>(blockIdx.y);
    if (origin >= patches)
        return;
    double sum = stayers[static_cast<size_t>(c) * patches + origin];
    // Neighboring threads read neighboring origins at each destination.
    for (int dest = 0; dest < patches; ++dest) {
        if (dest != origin) {
            const int t = origin - (origin > dest);
            sum += travelers[(static_cast<size_t>(dest) * compartments + c) * (patches - 1) + t];
        }
    }
    residents[static_cast<size_t>(c) * patches + origin] = sum;
}

template <int G>
cudaError_t launch_home(double* totals, const double* beta, const double* rate_e, const double* rate_i, int patches,
                        int steps, double dt, cudaStream_t stream)
{
    const int threads = patches <= 32 ? 32 : 64;
    seir_totals_allpatches_ct_kernel<G, false><<<(patches + threads - 1) / threads, threads, 0, stream>>>(
        totals, nullptr, beta, rate_e, rate_i, patches, steps, dt);
    return cudaGetLastError();
}
} // namespace

extern "C" cudaError_t launch_seir_stage_aligned_rk4_steps(double* totals, double* travelers, double* stage_lambda,
                                                           const double* contact_beta, const double* rate_exposed,
                                                           const double* rate_infected, int patches,
                                                           int travelers_per_patch, int groups, int steps, double dt,
                                                           cudaStream_t stream)
{
    if (totals == nullptr || travelers == nullptr || stage_lambda == nullptr || contact_beta == nullptr ||
        rate_exposed == nullptr || rate_infected == nullptr || patches <= 0 || travelers_per_patch <= 0 ||
        patches > 65535 || steps <= 0 || !std::isfinite(dt) || dt <= 0.0) {
        return cudaErrorInvalidValue;
    }

    switch (groups) {
    case 1:
        return launch_step<1>(totals, travelers, stage_lambda, contact_beta, rate_exposed, rate_infected, patches,
                              travelers_per_patch, steps, dt, stream);
    case 3:
        return launch_step<3>(totals, travelers, stage_lambda, contact_beta, rate_exposed, rate_infected, patches,
                              travelers_per_patch, steps, dt, stream);
    case 6:
        return launch_step<6>(totals, travelers, stage_lambda, contact_beta, rate_exposed, rate_infected, patches,
                              travelers_per_patch, steps, dt, stream);
    case 8:
        return launch_step<8>(totals, travelers, stage_lambda, contact_beta, rate_exposed, rate_infected, patches,
                              travelers_per_patch, steps, dt, stream);
    default:
        return cudaErrorInvalidValue;
    }
}

extern "C" cudaError_t launch_seir_daily_home(double* totals, const double* beta, const double* rate_e,
                                              const double* rate_i, int patches, int groups, int steps, double dt,
                                              cudaStream_t stream)
{
    if (!totals || !beta || !rate_e || !rate_i || patches < 2 || patches > 65535 || steps < 1 || !std::isfinite(dt) ||
        dt <= 0) {
        return cudaErrorInvalidValue;
    }
    switch (groups) {
    case 1:
        return launch_home<1>(totals, beta, rate_e, rate_i, patches, steps, dt, stream);
    case 3:
        return launch_home<3>(totals, beta, rate_e, rate_i, patches, steps, dt, stream);
    case 6:
        return launch_home<6>(totals, beta, rate_e, rate_i, patches, steps, dt, stream);
    case 8:
        return launch_home<8>(totals, beta, rate_e, rate_i, patches, steps, dt, stream);
    default:
        return cudaErrorInvalidValue;
    }
}

extern "C" cudaError_t launch_seir_daily_departure(double* totals, double* travelers, double* resident_snapshot,
                                                   const double* ht, int patches, int groups, cudaStream_t stream)
{
    if (!totals || !travelers || !resident_snapshot || !ht || patches < 2 || patches > 65535 ||
        (groups != 1 && groups != 3 && groups != 6 && groups != 8))
        return cudaErrorInvalidValue;
    auto status = cudaMemcpyAsync(resident_snapshot, totals, static_cast<size_t>(4 * groups) * patches * sizeof(double),
                                  cudaMemcpyDeviceToDevice, stream);
    if (status != cudaSuccess)
        return status;
    daily_departure_kernel<<<dim3(patches, 4 * groups), 256, 0, stream>>>(resident_snapshot, ht, totals, travelers,
                                                                          patches, 4 * groups);
    return cudaGetLastError();
}

extern "C" cudaError_t launch_seir_daily_return(double* totals, const double* travelers, double* stayers, int patches,
                                                int groups, cudaStream_t stream)
{
    if (!totals || !travelers || !stayers || patches < 2 || patches > 65535 ||
        (groups != 1 && groups != 3 && groups != 6 && groups != 8))
        return cudaErrorInvalidValue;
    daily_stayers_kernel<<<dim3(patches, 4 * groups), 256, 0, stream>>>(totals, travelers, stayers, patches,
                                                                        4 * groups);
    auto status = cudaGetLastError();
    if (status != cudaSuccess)
        return status;
    daily_return_kernel<<<dim3((patches + 255) / 256, 4 * groups), 256, 0, stream>>>(stayers, travelers, totals,
                                                                                     patches, 4 * groups);
    return cudaGetLastError();
}
