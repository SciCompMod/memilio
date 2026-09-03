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

#include <cstddef>

namespace
{

constexpr int block_size = 256;

template <int G>
__global__ void seir_totals_allpatches_ct_kernel(double* __restrict__ totals, double* __restrict__ stage_lambda,
                                                 const double* __restrict__ contact_beta,
                                                 const double* __restrict__ rate_exposed,
                                                 const double* __restrict__ rate_infected, int patches, double dt)
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
#pragma unroll
    for (int stage = 0; stage < 4; ++stage) {
#pragma unroll
        for (int target = 0; target < G; ++target) {
            lambda[target] = 0.0;
#pragma unroll
            for (int source = 0; source < G; ++source) {
                lambda[target] =
                    fma(parameters[target * G + source], i[source] * inverse_population[source], lambda[target]);
            }
            stage_lambda[(stage * G + target) * patches + patch] = lambda[target];
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
                                   int patches, int travelers_per_patch, double dt)
{
    __shared__ double parameters[6 * G];
    const int patch = static_cast<int>(blockIdx.y);
    for (int index = static_cast<int>(threadIdx.x); index < 4 * G; index += static_cast<int>(blockDim.x)) {
        parameters[index] = stage_lambda[index * patches + patch];
    }
    if (static_cast<int>(threadIdx.x) < G) {
        parameters[4 * G + threadIdx.x] = rate_exposed[threadIdx.x];
        parameters[5 * G + threadIdx.x] = rate_infected[threadIdx.x];
    }
    __syncthreads();

    const int traveler = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (traveler >= travelers_per_patch) {
        return;
    }

    const size_t patch_offset = static_cast<size_t>(patch) * 4 * G * travelers_per_patch;
#pragma unroll
    for (int group = 0; group < G; ++group) {
        const size_t susceptible_index = patch_offset + static_cast<size_t>(4 * group) * travelers_per_patch + traveler;
        const size_t exposed_index     = susceptible_index + travelers_per_patch;
        const size_t infected_index    = exposed_index + travelers_per_patch;
        const size_t recovered_index   = infected_index + travelers_per_patch;
        double susceptible             = travelers[susceptible_index];
        double exposed                 = travelers[exposed_index];
        double infected                = travelers[infected_index];
        double recovered               = travelers[recovered_index];
        const double lambda[4]         = {parameters[group], parameters[G + group], parameters[2 * G + group],
                                  parameters[3 * G + group]};
        integrate_traveler_group(susceptible, exposed, infected, recovered, lambda, parameters[4 * G + group],
                                 parameters[5 * G + group], dt);
        travelers[susceptible_index] = susceptible;
        travelers[exposed_index]     = exposed;
        travelers[infected_index]    = infected;
        travelers[recovered_index]   = recovered;
    }
}

template <int G>
cudaError_t launch_step(double* totals, double* travelers, double* stage_lambda, const double* contact_beta,
                        const double* rate_exposed, const double* rate_infected, int patches, int travelers_per_patch,
                        double dt, cudaStream_t stream)
{
    int totals_threads = 32;
    while (totals_threads < patches && totals_threads < block_size) {
        totals_threads *= 2;
    }
    const int totals_blocks = (patches + totals_threads - 1) / totals_threads;
    seir_totals_allpatches_ct_kernel<G><<<totals_blocks, totals_threads, 0, stream>>>(
        totals, stage_lambda, contact_beta, rate_exposed, rate_infected, patches, dt);
    auto status = cudaGetLastError();
    if (status != cudaSuccess) {
        return status;
    }

    int traveler_threads = 32;
    while (traveler_threads < travelers_per_patch && traveler_threads < block_size) {
        traveler_threads *= 2;
    }
    const dim3 traveler_grid(static_cast<unsigned int>((travelers_per_patch + traveler_threads - 1) / traveler_threads),
                             static_cast<unsigned int>(patches));
    seir_traveler_allpatches_ct_kernel<G><<<traveler_grid, traveler_threads, 0, stream>>>(
        travelers, stage_lambda, rate_exposed, rate_infected, patches, travelers_per_patch, dt);
    return cudaGetLastError();
}

} // namespace

extern "C" cudaError_t launch_seir_stage_aligned_rk4_step(double* totals, double* travelers, double* stage_lambda,
                                                          const double* contact_beta, const double* rate_exposed,
                                                          const double* rate_infected, int patches,
                                                          int travelers_per_patch, int groups, double dt,
                                                          cudaStream_t stream)
{
    if (totals == nullptr || travelers == nullptr || stage_lambda == nullptr || contact_beta == nullptr ||
        rate_exposed == nullptr || rate_infected == nullptr || patches <= 0 || travelers_per_patch <= 0 ||
        patches > 65535 || dt <= 0.0) {
        return cudaErrorInvalidValue;
    }

    switch (groups) {
    case 1:
        return launch_step<1>(totals, travelers, stage_lambda, contact_beta, rate_exposed, rate_infected, patches,
                              travelers_per_patch, dt, stream);
    case 3:
        return launch_step<3>(totals, travelers, stage_lambda, contact_beta, rate_exposed, rate_infected, patches,
                              travelers_per_patch, dt, stream);
    case 6:
        return launch_step<6>(totals, travelers, stage_lambda, contact_beta, rate_exposed, rate_infected, patches,
                              travelers_per_patch, dt, stream);
    case 8:
        return launch_step<8>(totals, travelers, stage_lambda, contact_beta, rate_exposed, rate_infected, patches,
                              travelers_per_patch, dt, stream);
    default:
        return cudaErrorInvalidValue;
    }
}
