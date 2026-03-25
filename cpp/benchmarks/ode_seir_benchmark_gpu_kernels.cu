/*
* Copyright (C) 2020-2025 MEmilio
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

/**
 *
 * Intentionally contains NO MEmilio headers!!!!
 * so nvcc never encounters the complex variadic-template code that breaks it.
 */

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/** Maximum supported compartments per commuter (NC = 4 × G, G <= 16). */
static constexpr int GPU_MAX_NC = 64;

/** Threads per CUDA block. */
static constexpr int GPU_BLOCK_SIZE = 256;

// ---------------------------------------------------------------------------
// Device helper: SEIR right-hand side (register-only, no global memory)
// ---------------------------------------------------------------------------

/**
 * @brief Evaluate the SEIR ODE RHS for one commuter.
 *
 * All data lives in registers; the function is forcibly inlined so that
 * register allocation happens inside the caller kernel.
 *
 * @param Xc     Current state vector [NC] (register array, read-only here).
 * @param k_out  Output derivative vector [NC] (register array).
 * @param lam    Force-of-infection per age group [G] (device pointer).
 * @param rE     E->I rate per age group [G] (device pointer).
 * @param rI     I->R rate per age group [G] (device pointer).
 * @param iS/E/I/R Flat compartment indices per age group [G] (device pointer).
 * @param G      Number of age groups.
 */
__device__ __forceinline__ void seir_rhs_dev(const double* __restrict__ Xc, double* __restrict__ k_out,
                                             const double* __restrict__ lam, const double* __restrict__ rE,
                                             const double* __restrict__ rI, const int* __restrict__ iS,
                                             const int* __restrict__ iE, const int* __restrict__ iI,
                                             const int* __restrict__ iR, int G)
{
    for (int g = 0; g < G; ++g) {
        const double fSE = lam[g] * Xc[iS[g]];
        const double fEI = rE[g] * Xc[iE[g]];
        const double fIR = rI[g] * Xc[iI[g]];
        k_out[iS[g]]     = -fSE;
        k_out[iE[g]]     = fSE - fEI;
        k_out[iI[g]]     = fEI - fIR;
        k_out[iR[g]]     = fIR;
    }
}

/**
 * @brief Compute force-of-infection lambda[i] from population y.
 *
 * lambda[i] = sum_j  contact[i*G+j] * beta[i] / N_j * I_j
 *
 * @param y        State vector [NC] (device pointer).
 * @param lam      Output lambda [G] (device pointer).
 * @param d_contact Row-major contact matrix [G*G].
 * @param d_beta   Transmission probability per age group [G].
 */
__device__ __forceinline__ void seir_compute_lambda_dev(const double* __restrict__ y, double* __restrict__ lam,
                                                        const double* __restrict__ d_contact,
                                                        const double* __restrict__ d_beta, const int* __restrict__ iS,
                                                        const int* __restrict__ iE, const int* __restrict__ iI,
                                                        const int* __restrict__ iR, int G)
{
    for (int i = 0; i < G; ++i)
        lam[i] = 0.0;
    for (int j = 0; j < G; ++j) {
        const double Nj    = y[iS[j]] + y[iE[j]] + y[iI[j]] + y[iR[j]];
        const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;
        const double Ij    = y[iI[j]];
        for (int i = 0; i < G; ++i)
            lam[i] += d_contact[i * G + j] * d_beta[i] * divNj * Ij;
    }
}

/**
 * @brief RK4 commuter update – one CUDA thread per commuter group.
 *
 * Register budget: 4 × NC doubles (Xc, k, acc, tmp).
 * For NC = 64 (G = 16): 256 doubles -> 2048 bytes of registers per thread.
 *
 * Memory access pattern: thread cg and thread cg+1 read the same compartment
 * at addresses [c * num_commuters + cg] and [c * num_commuters + cg+1]
 * -> consecutive -> fully coalesced.
 */
__global__ void seir_commuter_rk4_kernel(double* __restrict__ d_mobile,
                                         const double* __restrict__ d_lambda_stages, // [4 * G]
                                         const double* __restrict__ d_rE, // [G]
                                         const double* __restrict__ d_rI, // [G]
                                         const int* __restrict__ d_iS, // [G]
                                         const int* __restrict__ d_iE, // [G]
                                         const int* __restrict__ d_iI, // [G]
                                         const int* __restrict__ d_iR, // [G]
                                         int NC, int G, double dt, int num_commuters)
{
    const int cg = static_cast<int>(blockIdx.x) * GPU_BLOCK_SIZE + static_cast<int>(threadIdx.x);
    if (cg >= num_commuters)
        return;

    // ---- Load state into registers ----
    double Xc[GPU_MAX_NC];
    for (int c = 0; c < NC; ++c)
        Xc[c] = d_mobile[c * num_commuters + cg];

    // Running accumulator: k1 + 2*k2 + 2*k3  (k4 added at the end)
    double acc[GPU_MAX_NC];
    double k[GPU_MAX_NC];
    double tmp[GPU_MAX_NC];

    // ---- Stage 1: k1 ----
    seir_rhs_dev(Xc, k, d_lambda_stages + 0 * G, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);
    for (int c = 0; c < NC; ++c) {
        acc[c] = k[c]; // acc = k1
        tmp[c] = Xc[c] + (dt * 0.5) * k[c]; // tmp = Xc + (h/2)*k1
    }

    // ---- Stage 2: k2 ----
    seir_rhs_dev(tmp, k, d_lambda_stages + 1 * G, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);
    for (int c = 0; c < NC; ++c) {
        acc[c] += 2.0 * k[c]; // acc += 2*k2
        tmp[c] = Xc[c] + (dt * 0.5) * k[c]; // tmp = Xc + (h/2)*k2
    }

    // ---- Stage 3: k3 ----
    seir_rhs_dev(tmp, k, d_lambda_stages + 2 * G, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);
    for (int c = 0; c < NC; ++c) {
        acc[c] += 2.0 * k[c]; // acc += 2*k3
        tmp[c] = Xc[c] + dt * k[c]; // tmp = Xc + h*k3
    }

    // ---- Stage 4: k4 ----
    seir_rhs_dev(tmp, k, d_lambda_stages + 3 * G, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // ---- Write back ----
    const double c6 = dt / 6.0;
    for (int c = 0; c < NC; ++c)
        d_mobile[c * num_commuters + cg] = Xc[c] + c6 * (acc[c] + k[c]);
}

/**
 * @brief Euler commuter update – one CUDA thread per commuter group.
 *
 * Register budget: 1 × NC doubles.  Fluxes are computed from the old state
 * before any compartment in the same age group is updated.
 */
__global__ void seir_commuter_euler_kernel(double* __restrict__ d_mobile,
                                           const double* __restrict__ d_lambda, // [G]
                                           const double* __restrict__ d_rE, const double* __restrict__ d_rI,
                                           const int* __restrict__ d_iS, const int* __restrict__ d_iE,
                                           const int* __restrict__ d_iI, const int* __restrict__ d_iR, int NC, int G,
                                           double dt, int num_commuters)
{
    const int cg = static_cast<int>(blockIdx.x) * GPU_BLOCK_SIZE + static_cast<int>(threadIdx.x);
    if (cg >= num_commuters)
        return;

    double Xc[GPU_MAX_NC];
    for (int c = 0; c < NC; ++c)
        Xc[c] = d_mobile[c * num_commuters + cg];

    for (int g = 0; g < G; ++g) {
        const double fSE = d_lambda[g] * Xc[d_iS[g]];
        const double fEI = d_rE[g] * Xc[d_iE[g]];
        const double fIR = d_rI[g] * Xc[d_iI[g]];
        Xc[d_iS[g]] += dt * (-fSE);
        Xc[d_iE[g]] += dt * (fSE - fEI);
        Xc[d_iI[g]] += dt * (fEI - fIR);
        Xc[d_iR[g]] += dt * fIR;
    }

    for (int c = 0; c < NC; ++c)
        d_mobile[c * num_commuters + cg] = Xc[c];
}

/**
 * @brief RK4 for aggregate totals + lambda stages – runs on a single GPU thread.
 *
 * This kernel eliminates the per-step H2D lambda transfer and CPU-side
 * derivative computation from the stage-aligned benchmark.
 *
 * Inputs/Outputs:
 *   d_z              [NC]   in-out: totals state vector, advanced by one step.
 *   d_lambda_stages  [4*G]  output: lambda[stage*G .. stage*G+G-1] per RK4 stage.
 *   d_contact        [G*G]  row-major contact matrix (constant across steps).
 *   d_beta           [G]    transmission probability (constant).
 *
 * Register budget: 7×NC doubles (z, k1..k4, tmp) + 1×G lambda 
 * For NC=64, G=16: 448 + 16 = 464 doubles 
 */
__global__ void seir_totals_rk4_lambda_kernel(double* __restrict__ d_z, double* __restrict__ d_lambda_stages,
                                              const double* __restrict__ d_contact, const double* __restrict__ d_beta,
                                              const double* __restrict__ d_rE, const double* __restrict__ d_rI,
                                              const int* __restrict__ d_iS, const int* __restrict__ d_iE,
                                              const int* __restrict__ d_iI, const int* __restrict__ d_iR, int NC, int G,
                                              double dt)
{
    if (threadIdx.x != 0 || blockIdx.x != 0)
        return;

    double z[GPU_MAX_NC];
    double k1[GPU_MAX_NC];
    double k2[GPU_MAX_NC];
    double k3[GPU_MAX_NC];
    double k4[GPU_MAX_NC];
    double tmp[GPU_MAX_NC];
    double lam[GPU_MAX_NC / 4]; // G <= 16

    for (int c = 0; c < NC; ++c)
        z[c] = d_z[c];

    // ---- Stage 1 ----
    seir_compute_lambda_dev(z, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR, G);
    for (int g = 0; g < G; ++g)
        d_lambda_stages[0 * G + g] = lam[g];
    seir_rhs_dev(z, k1, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // ---- Stage 2 ----
    for (int c = 0; c < NC; ++c)
        tmp[c] = z[c] + 0.5 * dt * k1[c];
    seir_compute_lambda_dev(tmp, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR, G);
    for (int g = 0; g < G; ++g)
        d_lambda_stages[1 * G + g] = lam[g];
    seir_rhs_dev(tmp, k2, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // ---- Stage 3 ----
    for (int c = 0; c < NC; ++c)
        tmp[c] = z[c] + 0.5 * dt * k2[c];
    seir_compute_lambda_dev(tmp, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR, G);
    for (int g = 0; g < G; ++g)
        d_lambda_stages[2 * G + g] = lam[g];
    seir_rhs_dev(tmp, k3, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // ---- Stage 4 ----
    for (int c = 0; c < NC; ++c)
        tmp[c] = z[c] + dt * k3[c];
    seir_compute_lambda_dev(tmp, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR, G);
    for (int g = 0; g < G; ++g)
        d_lambda_stages[3 * G + g] = lam[g];
    seir_rhs_dev(tmp, k4, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // ---- RK4 update ----
    const double c6 = dt / 6.0;
    for (int c = 0; c < NC; ++c)
        d_z[c] = z[c] + c6 * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]);
}

/** Compile-time SEIR RHS. */
template <int G_>
__device__ __forceinline__ void
seir_rhs_ct(const double* __restrict__ Xc, double* __restrict__ k_out, const double* __restrict__ lam,
            const double* __restrict__ rE, const double* __restrict__ rI, const int* __restrict__ iS,
            const int* __restrict__ iE, const int* __restrict__ iI, const int* __restrict__ iR)
{
#pragma unroll
    for (int g = 0; g < G_; ++g) {
        const double fSE = lam[g] * Xc[iS[g]];
        const double fEI = rE[g] * Xc[iE[g]];
        const double fIR = rI[g] * Xc[iI[g]];
        k_out[iS[g]]     = -fSE;
        k_out[iE[g]]     = fSE - fEI;
        k_out[iI[g]]     = fEI - fIR;
        k_out[iR[g]]     = fIR;
    }
}

/** Compile-time lambda computation (G×G contact loop unrolled). */
template <int G_>
__device__ __forceinline__ void seir_lambda_ct(const double* __restrict__ y, double* __restrict__ lam,
                                               const double* __restrict__ d_contact, const double* __restrict__ d_beta,
                                               const int* __restrict__ iS, const int* __restrict__ iE,
                                               const int* __restrict__ iI, const int* __restrict__ iR)
{
#pragma unroll
    for (int i = 0; i < G_; ++i)
        lam[i] = 0.0;
#pragma unroll
    for (int j = 0; j < G_; ++j) {
        const double Nj    = y[iS[j]] + y[iE[j]] + y[iI[j]] + y[iR[j]];
        const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;
        const double Ij    = y[iI[j]];
#pragma unroll
        for (int i = 0; i < G_; ++i)
            lam[i] += d_contact[i * G_ + j] * d_beta[i] * divNj * Ij;
    }
}

/**
 * @brief Compile-time totals RK4 + lambda for all patches.
 *
 * G_ is a compile-time constant -> all local arrays are register-only,
 * all loops are unrolled.  NC_ = 4*G_ doubles per thread in registers.
 * The P and dt parameters remain runtime values (block count is computed
 * by the host wrapper).
 */
template <int G_>
__global__ void seir_totals_allpatches_ct_kernel(double* __restrict__ d_z, double* __restrict__ d_lambda_stages,
                                                 const double* __restrict__ d_contact,
                                                 const double* __restrict__ d_beta, const double* __restrict__ d_rE,
                                                 const double* __restrict__ d_rI, const int* __restrict__ d_iS,
                                                 const int* __restrict__ d_iE, const int* __restrict__ d_iI,
                                                 const int* __restrict__ d_iR, int P, double dt)
{
    constexpr int NC_ = 4 * G_;
    const int p       = static_cast<int>(blockIdx.x) * GPU_BLOCK_SIZE + static_cast<int>(threadIdx.x);
    if (p >= P)
        return;
    double z[NC_]; // original state, held as base for all stage inputs
    double acc[NC_]; // running accumulator: z0 + dt*(k1 + 2k2 + 2k3 + k4)/6
    double k[NC_]; // current-stage derivative (reused)
    double tmp[NC_]; // next-stage input (reused)
    double lam[G_]; // lambda (reused per stage)

#pragma unroll
    for (int c = 0; c < NC_; ++c)
        z[c] = acc[c] = d_z[c * P + p];

    // Stage 1: k1, lambda_0
    seir_lambda_ct<G_>(z, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int g = 0; g < G_; ++g)
        d_lambda_stages[(0 * G_ + g) * P + p] = lam[g];
    seir_rhs_ct<G_>(z, k, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int c = 0; c < NC_; ++c) {
        acc[c] += (dt / 6.0) * k[c]; // acc += dt/6 * k1
        tmp[c] = z[c] + 0.5 * dt * k[c]; // input for stage 2
    }

    // Stage 2: k2, lambda_1
    seir_lambda_ct<G_>(tmp, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int g = 0; g < G_; ++g)
        d_lambda_stages[(1 * G_ + g) * P + p] = lam[g];
    seir_rhs_ct<G_>(tmp, k, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int c = 0; c < NC_; ++c) {
        acc[c] += (dt / 3.0) * k[c]; // acc += dt/3 * k2
        tmp[c] = z[c] + 0.5 * dt * k[c]; // input for stage 3
    }

    // Stage 3: k3, lambda_2
    seir_lambda_ct<G_>(tmp, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int g = 0; g < G_; ++g)
        d_lambda_stages[(2 * G_ + g) * P + p] = lam[g];
    seir_rhs_ct<G_>(tmp, k, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int c = 0; c < NC_; ++c) {
        acc[c] += (dt / 3.0) * k[c]; // acc += dt/3 * k3
        tmp[c] = z[c] + dt * k[c]; // input for stage 4
    }

    // Stage 4: k4, lambda_3
    seir_lambda_ct<G_>(tmp, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int g = 0; g < G_; ++g)
        d_lambda_stages[(3 * G_ + g) * P + p] = lam[g];
    seir_rhs_ct<G_>(tmp, k, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int c = 0; c < NC_; ++c)
        acc[c] += (dt / 6.0) * k[c]; // acc += dt/6 * k4  -> acc = z0 + dt*(k1+2k2+2k3+k4)/6

#pragma unroll
    for (int c = 0; c < NC_; ++c)
        d_z[c * P + p] = acc[c];
}

/**
 * @brief Compile-time commuter RK4 for all (patch, commuter) pairs.
 */
template <int G_>
__global__ void
seir_commuter_allpatches_ct_kernel(double* __restrict__ d_mobile, const double* __restrict__ d_lambda_stages,
                                   const double* __restrict__ d_rE, const double* __restrict__ d_rI,
                                   const int* __restrict__ d_iS, const int* __restrict__ d_iE,
                                   const int* __restrict__ d_iI, const int* __restrict__ d_iR, int P, int N, double dt)
{
    constexpr int NC_ = 4 * G_;
    const int p       = static_cast<int>(blockIdx.y);
    const int cg      = static_cast<int>(blockIdx.x) * GPU_BLOCK_SIZE + static_cast<int>(threadIdx.x);
    if (cg >= N || p >= P)
        return;

    // Accumulator pattern: 4 arrays instead of 7.
    // Xc[NC_] + acc[NC_] + k[NC_] + tmp[NC_] + lam[G_] = 17G doubles.
    double Xc[NC_]; // initial state (base for all stage inputs)
    double acc[NC_]; // accumulator: starts = Xc, += weighted k per stage
    double k[NC_]; // current derivative (reused across stages)
    double tmp[NC_]; // stage input (reused across stages)
    double lam[G_]; // lambda (reused per stage)

#pragma unroll
    for (int c = 0; c < NC_; ++c)
        Xc[c] = acc[c] = d_mobile[(p * NC_ + c) * N + cg];

// Stage 1
#pragma unroll
    for (int g = 0; g < G_; ++g)
        lam[g] = d_lambda_stages[(0 * G_ + g) * P + p];
    seir_rhs_ct<G_>(Xc, k, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int c = 0; c < NC_; ++c) {
        acc[c] += (dt / 6.0) * k[c];
        tmp[c] = Xc[c] + 0.5 * dt * k[c];
    }

// Stage 2
#pragma unroll
    for (int g = 0; g < G_; ++g)
        lam[g] = d_lambda_stages[(1 * G_ + g) * P + p];
    seir_rhs_ct<G_>(tmp, k, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int c = 0; c < NC_; ++c) {
        acc[c] += (dt / 3.0) * k[c];
        tmp[c] = Xc[c] + 0.5 * dt * k[c];
    }

// Stage 3
#pragma unroll
    for (int g = 0; g < G_; ++g)
        lam[g] = d_lambda_stages[(2 * G_ + g) * P + p];
    seir_rhs_ct<G_>(tmp, k, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int c = 0; c < NC_; ++c) {
        acc[c] += (dt / 3.0) * k[c];
        tmp[c] = Xc[c] + dt * k[c];
    }

// Stage 4
#pragma unroll
    for (int g = 0; g < G_; ++g)
        lam[g] = d_lambda_stages[(3 * G_ + g) * P + p];
    seir_rhs_ct<G_>(tmp, k, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR);
#pragma unroll
    for (int c = 0; c < NC_; ++c)
        acc[c] += (dt / 6.0) * k[c]; // acc = Xc + dt*(k1+2k2+2k3+k4)/6

#pragma unroll
    for (int c = 0; c < NC_; ++c)
        d_mobile[(p * NC_ + c) * N + cg] = acc[c];
}

/**
 * @brief RK4 + lambda for all patches simultaneously – one thread per patch.
 *
 * Launched with <<<ceil(P/BLOCK), BLOCK>>>.
 * Each thread owns one column of d_z (a single patch's NC-vector) and
 * performs the full 4-stage RK4 integration, writing lambda for every
 * stage into d_lambda_stages so the commuter kernel can
 * read it without any CPU involvement.
 */
__global__ void
seir_totals_rk4_lambda_allpatches_kernel(double* __restrict__ d_z, // [NC × P], z[c*P + p]
                                         double* __restrict__ d_lambda_stages, // [4G × P], lam[(s*G+g)*P + p]
                                         const double* __restrict__ d_contact, // [G×G] shared
                                         const double* __restrict__ d_beta, const double* __restrict__ d_rE,
                                         const double* __restrict__ d_rI, const int* __restrict__ d_iS,
                                         const int* __restrict__ d_iE, const int* __restrict__ d_iI,
                                         const int* __restrict__ d_iR, int NC, int G, int P, double dt)
{
    const int p = static_cast<int>(blockIdx.x) * GPU_BLOCK_SIZE + static_cast<int>(threadIdx.x);
    if (p >= P)
        return;

    double z[GPU_MAX_NC];
    double k1[GPU_MAX_NC];
    double k2[GPU_MAX_NC];
    double k3[GPU_MAX_NC];
    double k4[GPU_MAX_NC];
    double tmp[GPU_MAX_NC];
    double lam[GPU_MAX_NC / 4]; // G <= 16

    for (int c = 0; c < NC; ++c)
        z[c] = d_z[c * P + p];

    // --- Stage 1 ---
    seir_compute_lambda_dev(z, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR, G);
    for (int g = 0; g < G; ++g)
        d_lambda_stages[(0 * G + g) * P + p] = lam[g];
    seir_rhs_dev(z, k1, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // --- Stage 2 ---
    for (int c = 0; c < NC; ++c)
        tmp[c] = z[c] + 0.5 * dt * k1[c];
    seir_compute_lambda_dev(tmp, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR, G);
    for (int g = 0; g < G; ++g)
        d_lambda_stages[(1 * G + g) * P + p] = lam[g];
    seir_rhs_dev(tmp, k2, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // --- Stage 3 ---
    for (int c = 0; c < NC; ++c)
        tmp[c] = z[c] + 0.5 * dt * k2[c];
    seir_compute_lambda_dev(tmp, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR, G);
    for (int g = 0; g < G; ++g)
        d_lambda_stages[(2 * G + g) * P + p] = lam[g];
    seir_rhs_dev(tmp, k3, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // --- Stage 4 ---
    for (int c = 0; c < NC; ++c)
        tmp[c] = z[c] + dt * k3[c];
    seir_compute_lambda_dev(tmp, lam, d_contact, d_beta, d_iS, d_iE, d_iI, d_iR, G);
    for (int g = 0; g < G; ++g)
        d_lambda_stages[(3 * G + g) * P + p] = lam[g];
    seir_rhs_dev(tmp, k4, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // --- RK4 update ---
    const double c6 = dt / 6.0;
    for (int c = 0; c < NC; ++c)
        d_z[c * P + p] = z[c] + c6 * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]);
}

/**
 * @brief RK4 commuter update for all patches simultaneously.
 *
 * 2-D grid:  blockIdx.y = patch index (0..P-1)
 *            blockIdx.x * BLOCK + threadIdx.x = commuter index (0..N-1)
 *
 * Each thread handles exactly one (patch, commuter) pair and performs the
 * full 4-stage RK4 step, reading lambda from d_lambda_stages which was
 * written by seir_totals_rk4_lambda_allpatches_kernel on the same stream.
 */
__global__ void seir_commuter_rk4_allpatches_kernel(double* __restrict__ d_mobile, // [(c*P + p)*N + cg]
                                                    const double* __restrict__ d_lambda_stages, // [(s*G + g)*P + p]
                                                    const double* __restrict__ d_rE, const double* __restrict__ d_rI,
                                                    const int* __restrict__ d_iS, const int* __restrict__ d_iE,
                                                    const int* __restrict__ d_iI, const int* __restrict__ d_iR, int NC,
                                                    int G, int P, int N, double dt)
{
    const int p  = static_cast<int>(blockIdx.y);
    const int cg = static_cast<int>(blockIdx.x) * GPU_BLOCK_SIZE + static_cast<int>(threadIdx.x);
    if (cg >= N || p >= P)
        return;

    double Xc[GPU_MAX_NC];
    double k1[GPU_MAX_NC];
    double k2[GPU_MAX_NC];
    double k3[GPU_MAX_NC];
    double k4[GPU_MAX_NC];
    double tmp[GPU_MAX_NC];
    double lam[GPU_MAX_NC / 4]; // G <= 16

    for (int c = 0; c < NC; ++c)
        Xc[c] = d_mobile[(p * NC + c) * N + cg];

    // --- Stage 1 ---
    for (int g = 0; g < G; ++g)
        lam[g] = d_lambda_stages[(0 * G + g) * P + p];
    seir_rhs_dev(Xc, k1, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // --- Stage 2 ---
    for (int c = 0; c < NC; ++c)
        tmp[c] = Xc[c] + 0.5 * dt * k1[c];
    for (int g = 0; g < G; ++g)
        lam[g] = d_lambda_stages[(1 * G + g) * P + p];
    seir_rhs_dev(tmp, k2, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // --- Stage 3 ---
    for (int c = 0; c < NC; ++c)
        tmp[c] = Xc[c] + 0.5 * dt * k2[c];
    for (int g = 0; g < G; ++g)
        lam[g] = d_lambda_stages[(2 * G + g) * P + p];
    seir_rhs_dev(tmp, k3, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // --- Stage 4 ---
    for (int c = 0; c < NC; ++c)
        tmp[c] = Xc[c] + dt * k3[c];
    for (int g = 0; g < G; ++g)
        lam[g] = d_lambda_stages[(3 * G + g) * P + p];
    seir_rhs_dev(tmp, k4, lam, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, G);

    // --- RK4 update ---
    const double c6 = dt / 6.0;
    for (int c = 0; c < NC; ++c)
        d_mobile[(p * NC + c) * N + cg] = Xc[c] + c6 * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]);
}

// ---------------------------------------------------------------------------
// Kernel: Phi block-diagonal reconstruction  Xt = Phi * X0
// ---------------------------------------------------------------------------

/**
 * @brief Apply block-diagonal Phi reconstruction for all commuter groups.
 *
 * One CUDA thread per commuter group.  Each thread:
 *   1. Loads the full NC-component X0 vector into registers.
 *   2. Applies G independent 4×4 matrix-vector products (one per age group).
 *   3. Writes the result Xt back.
 */
__global__ void seir_phi_reconstruction_kernel(double* __restrict__ d_Xt, // [NC × num_commuters] col-major output
                                               const double* __restrict__ d_X0, // [NC × num_commuters] col-major input
                                               const double* __restrict__ d_Phi, // [G × 16] col-major 4×4 Phi blocks
                                               const int* __restrict__ d_iS, // [G]
                                               const int* __restrict__ d_iE, // [G]
                                               const int* __restrict__ d_iI, // [G]
                                               const int* __restrict__ d_iR, // [G]
                                               int NC, int G, int num_commuters)
{
    const int cg = static_cast<int>(blockIdx.x) * GPU_BLOCK_SIZE + static_cast<int>(threadIdx.x);
    if (cg >= num_commuters)
        return;

    // ---- Load X0 into registers (coalesced) ----
    double X0[GPU_MAX_NC];
    for (int c = 0; c < NC; ++c)
        X0[c] = d_X0[c * num_commuters + cg];

    // ---- Output in registers (zero-init) ----
    double Xt[GPU_MAX_NC];
    for (int c = 0; c < NC; ++c)
        Xt[c] = 0.0;

    // ---- Apply each 4×4 Phi block (all in registers / constant cache) ----
    for (int g = 0; g < G; ++g) {
        const double* P = d_Phi + g * 16; // col-major: P[row + col*4]
        const int iS = d_iS[g], iE = d_iE[g], iI = d_iI[g], iR = d_iR[g];

        const double x0S = X0[iS], x0E = X0[iE], x0I = X0[iI], x0R = X0[iR];

        // 4×4 col-major matrix-vector:  Xt_g = P_g * X0_g
        Xt[iS] = P[0] * x0S + P[4] * x0E + P[8] * x0I + P[12] * x0R;
        Xt[iE] = P[1] * x0S + P[5] * x0E + P[9] * x0I + P[13] * x0R;
        Xt[iI] = P[2] * x0S + P[6] * x0E + P[10] * x0I + P[14] * x0R;
        Xt[iR] = P[3] * x0S + P[7] * x0E + P[11] * x0I + P[15] * x0R;
    }

    // ---- Write Xt back (coalesced) ----
    for (int c = 0; c < NC; ++c)
        d_Xt[c * num_commuters + cg] = Xt[c];
}

// ---------------------------------------------------------------------------
// Public extern "C" wrappers (called from the host-side .cpp benchmark)
// ---------------------------------------------------------------------------

extern "C" {

/**
 * @brief Launch the RK4 commuter-update kernel.
 *
 * All pointers must already reside on the device.
 * The kernel is launched asynchronously on the default stream.
 */
void launch_seir_commuter_rk4(double* d_mobile, const double* d_lambda_stages, const double* d_rE, const double* d_rI,
                              const int* d_iS, const int* d_iE, const int* d_iI, const int* d_iR, int NC, int G,
                              double dt, int num_commuters)
{
    const int blocks = (num_commuters + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;
    seir_commuter_rk4_kernel<<<blocks, GPU_BLOCK_SIZE>>>(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR,
                                                         NC, G, dt, num_commuters);
}

/**
 * @brief Launch the Euler commuter-update kernel.
 *
 * All pointers must already reside on the device.
 * The kernel is launched asynchronously on the default stream.
 */
void launch_seir_commuter_euler(double* d_mobile, const double* d_lambda, const double* d_rE, const double* d_rI,
                                const int* d_iS, const int* d_iE, const int* d_iI, const int* d_iR, int NC, int G,
                                double dt, int num_commuters)
{
    const int blocks = (num_commuters + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;
    seir_commuter_euler_kernel<<<blocks, GPU_BLOCK_SIZE>>>(d_mobile, d_lambda, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, NC,
                                                           G, dt, num_commuters);
}

/**
 * @brief Launch the Phi block-diagonal reconstruction kernel.
 *
 * All pointers must already reside on the device.
 * The kernel is launched asynchronously on the default stream.
 */
void launch_seir_phi_reconstruction(double* d_Xt, const double* d_X0, const double* d_Phi, const int* d_iS,
                                    const int* d_iE, const int* d_iI, const int* d_iR, int NC, int G, int num_commuters)
{
    const int blocks = (num_commuters + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;
    seir_phi_reconstruction_kernel<<<blocks, GPU_BLOCK_SIZE>>>(d_Xt, d_X0, d_Phi, d_iS, d_iE, d_iI, d_iR, NC, G,
                                                               num_commuters);
}

/**
 * @brief Launch the totals-RK4 + lambda kernel (1 thread, async).
 *
 * Updates d_z by one RK4 step and writes 4×G lambda values to d_lambda_stages.
 * Always launched on the default stream, so the subsequent commuter kernel
 * automatically waits for it.
 */
void launch_seir_totals_rk4_lambda(double* d_z, double* d_lambda_stages, const double* d_contact, const double* d_beta,
                                   const double* d_rE, const double* d_rI, const int* d_iS, const int* d_iE,
                                   const int* d_iI, const int* d_iR, int NC, int G, double dt)
{
    seir_totals_rk4_lambda_kernel<<<1, 1>>>(d_z, d_lambda_stages, d_contact, d_beta, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR,
                                            NC, G, dt);
}

/**
 * @brief Launch totals RK4 + lambda kernel for all P patches in parallel.
 *
 * Launched with <<<ceil(P/BLOCK), BLOCK>>> on the default stream.
 * The subsequent commuter kernel must be launched on the same stream to
 * ensure it sees the freshly written d_lambda_stages.
 */
void launch_seir_totals_rk4_lambda_allpatches(double* d_z, double* d_lambda_stages, const double* d_contact,
                                              const double* d_beta, const double* d_rE, const double* d_rI,
                                              const int* d_iS, const int* d_iE, const int* d_iI, const int* d_iR,
                                              int NC, int G, int P, double dt)
{
    const int blocks = (P + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;
    // Dispatch to the compile-time-G kernel so all inner loops are fully
    // unrolled and the NC_-sized local arrays live in hardware registers.
    switch (G) {
    case 1:
        seir_totals_allpatches_ct_kernel<1><<<blocks, GPU_BLOCK_SIZE>>>(d_z, d_lambda_stages, d_contact, d_beta, d_rE,
                                                                        d_rI, d_iS, d_iE, d_iI, d_iR, P, dt);
        break;
    case 2:
        seir_totals_allpatches_ct_kernel<2><<<blocks, GPU_BLOCK_SIZE>>>(d_z, d_lambda_stages, d_contact, d_beta, d_rE,
                                                                        d_rI, d_iS, d_iE, d_iI, d_iR, P, dt);
        break;
    case 3:
        seir_totals_allpatches_ct_kernel<3><<<blocks, GPU_BLOCK_SIZE>>>(d_z, d_lambda_stages, d_contact, d_beta, d_rE,
                                                                        d_rI, d_iS, d_iE, d_iI, d_iR, P, dt);
        break;
    case 4:
        seir_totals_allpatches_ct_kernel<4><<<blocks, GPU_BLOCK_SIZE>>>(d_z, d_lambda_stages, d_contact, d_beta, d_rE,
                                                                        d_rI, d_iS, d_iE, d_iI, d_iR, P, dt);
        break;
    case 5:
        seir_totals_allpatches_ct_kernel<5><<<blocks, GPU_BLOCK_SIZE>>>(d_z, d_lambda_stages, d_contact, d_beta, d_rE,
                                                                        d_rI, d_iS, d_iE, d_iI, d_iR, P, dt);
        break;
    case 6:
        seir_totals_allpatches_ct_kernel<6><<<blocks, GPU_BLOCK_SIZE>>>(d_z, d_lambda_stages, d_contact, d_beta, d_rE,
                                                                        d_rI, d_iS, d_iE, d_iI, d_iR, P, dt);
        break;
    default:
        printf("launch_seir_totals_allpatches: G=%d not compiled in switch\n", G);
    }
}

/**
 * @brief Launch commuter RK4 kernel for all (patch, commuter) pairs in parallel.
 *
 * 2-D grid: dim3((ceil(N/BLOCK), P)  – one thread per (commuter, patch) pair.
 * Must be launched after launch_seir_totals_rk4_lambda_allpatches on the
 * same stream (d_lambda_stages dependency).
 */
void launch_seir_commuter_rk4_allpatches(double* d_mobile, const double* d_lambda_stages, const double* d_rE,
                                         const double* d_rI, const int* d_iS, const int* d_iE, const int* d_iI,
                                         const int* d_iR, int NC, int G, int P, int N, double dt)
{
    const int bx = (N + GPU_BLOCK_SIZE - 1) / GPU_BLOCK_SIZE;
    dim3 grid(static_cast<unsigned>(bx), static_cast<unsigned>(P));
    switch (G) {
    case 1:
        seir_commuter_allpatches_ct_kernel<1>
            <<<grid, GPU_BLOCK_SIZE>>>(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, P, N, dt);
        break;
    case 2:
        seir_commuter_allpatches_ct_kernel<2>
            <<<grid, GPU_BLOCK_SIZE>>>(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, P, N, dt);
        break;
    case 3:
        seir_commuter_allpatches_ct_kernel<3>
            <<<grid, GPU_BLOCK_SIZE>>>(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, P, N, dt);
        break;
    case 4:
        seir_commuter_allpatches_ct_kernel<4>
            <<<grid, GPU_BLOCK_SIZE>>>(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, P, N, dt);
        break;
    case 5:
        seir_commuter_allpatches_ct_kernel<5>
            <<<grid, GPU_BLOCK_SIZE>>>(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, P, N, dt);
        break;
    case 6:
        seir_commuter_allpatches_ct_kernel<6>
            <<<grid, GPU_BLOCK_SIZE>>>(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, P, N, dt);
        break;
    default:
        printf("launch_seir_commuter_allpatches: G=%d not compiled in switch\n", G);
    }
}

} // extern "C"
