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
 * This file is compiled by g++ (no nvcc).  All MEmilio template code is here.
 * The CUDA kernels are in ode_seir_benchmark_stage_aligned_gpu_kernels.cu.
 */

#include "benchmark/benchmark.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"
#include "ode_seir/model.h"

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------------------
// CUDA error-checking macro
// ---------------------------------------------------------------------------
#define CUDA_CHECK(call)                                                                                               \
    do {                                                                                                               \
        cudaError_t _e_ = (call);                                                                                      \
        if (_e_ != cudaSuccess) {                                                                                      \
            fprintf(stderr, "CUDA error at %s:%d – %s\n", __FILE__, __LINE__, cudaGetErrorString(_e_));                \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

// ---------------------------------------------------------------------------
// Declarations of the extern "C" kernel launchers (defined in the .cu file)
// ---------------------------------------------------------------------------
extern "C" {

void launch_seir_commuter_rk4(double* d_mobile, const double* d_lambda_stages, const double* d_rE, const double* d_rI,
                              const int* d_iS, const int* d_iE, const int* d_iI, const int* d_iR, int NC, int G,
                              double dt, int num_commuters);

void launch_seir_commuter_euler(double* d_mobile, const double* d_lambda, const double* d_rE, const double* d_rI,
                                const int* d_iS, const int* d_iE, const int* d_iI, const int* d_iR, int NC, int G,
                                double dt, int num_commuters);

void launch_seir_phi_reconstruction(double* d_Xt, const double* d_X0, const double* d_Phi, const int* d_iS,
                                    const int* d_iE, const int* d_iI, const int* d_iR, int NC, int G,
                                    int num_commuters);

void launch_seir_totals_rk4_lambda(double* d_z, double* d_lambda_stages, const double* d_contact, const double* d_beta,
                                   const double* d_rE, const double* d_rI, const int* d_iS, const int* d_iE,
                                   const int* d_iI, const int* d_iR, int NC, int G, double dt);

// All-patches variants (outer patch loop on GPU, 2-D grid for commuters)
void launch_seir_totals_rk4_lambda_allpatches(double* d_z, double* d_lambda_stages, const double* d_contact,
                                              const double* d_beta, const double* d_rE, const double* d_rI,
                                              const int* d_iS, const int* d_iE, const int* d_iI, const int* d_iR,
                                              int NC, int G, int P, double dt);

void launch_seir_commuter_rk4_allpatches(double* d_mobile, const double* d_lambda_stages, const double* d_rE,
                                         const double* d_rI, const int* d_iS, const int* d_iE, const int* d_iI,
                                         const int* d_iR, int NC, int G, int P, int N, double dt);

} // extern "C"

// ---------------------------------------------------------------------------
// Host helpers
// ---------------------------------------------------------------------------

namespace mio
{
namespace oseir
{

/** Cache of flat population indices for S, E, I, R per age group. */
struct SEIRIndexCache {
    std::vector<int> S, E, I, R;
    size_t G = 0;
};

inline SEIRIndexCache make_seir_index_cache(const mio::oseir::Model<ScalarType>& model)
{
    SEIRIndexCache ic;
    ic.G = static_cast<size_t>(model.parameters.get_num_groups());
    ic.S.resize(ic.G);
    ic.E.resize(ic.G);
    ic.I.resize(ic.G);
    ic.R.resize(ic.G);
    using IS = mio::oseir::InfectionState;
    for (size_t g = 0; g < ic.G; ++g) {
        auto Gg = mio::AgeGroup(static_cast<int>(g));
        ic.S[g] = model.populations.get_flat_index({Gg, IS::Susceptible});
        ic.E[g] = model.populations.get_flat_index({Gg, IS::Exposed});
        ic.I[g] = model.populations.get_flat_index({Gg, IS::Infected});
        ic.R[g] = model.populations.get_flat_index({Gg, IS::Recovered});
    }
    return ic;
}

} // namespace oseir
} // namespace mio

namespace mio
{
namespace benchmark_mio
{

using ModelType = mio::oseir::Model<ScalarType>;

const ScalarType t0_gpu    = 0.0;
const ScalarType t_max_gpu = 1.0;
const ScalarType dt_gpu    = 0.5;

const std::vector<int> commuter_group_counts_gpu = {256,  512,   1024, 2048, 4096,
                                                    8192, 16384, 32768}; // , 65536, 131072, 262144, 524288, 1048576};
const std::vector<int> age_group_counts_gpu      = {1, 2, 3, 4, 5, 6, 8};

/** Identical setup to the CPU reference benchmark. */
inline void setup_model_gpu(ModelType& model, size_t num_agegroups)
{
    const double total_population = 10000.0;
    for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, mio::oseir::InfectionState::Exposed}]   = 100.0 / static_cast<double>(num_agegroups);
        model.populations[{i, mio::oseir::InfectionState::Infected}]  = 100.0 / static_cast<double>(num_agegroups);
        model.populations[{i, mio::oseir::InfectionState::Recovered}] = 100.0 / static_cast<double>(num_agegroups);
        model.populations[{i, mio::oseir::InfectionState::Susceptible}] =
            total_population / static_cast<double>(num_agegroups) -
            model.populations[{i, mio::oseir::InfectionState::Exposed}] -
            model.populations[{i, mio::oseir::InfectionState::Infected}] -
            model.populations[{i, mio::oseir::InfectionState::Recovered}];
    }
    model.parameters.template set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model.parameters.template set<mio::oseir::TimeInfected<ScalarType>>(6.0);
    model.parameters.template set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);
    mio::ContactMatrixGroup& cm = model.parameters.template get<mio::oseir::ContactPatterns<ScalarType>>();
    cm[0].get_baseline().setConstant(2.7);
    model.apply_constraints();
}

/**
 * @brief GPU stage-aligned RK4 benchmark.
 *
 * Measured time: per-step totals ODE (CPU) + lambda computation (CPU) +
 *                tiny H2D upload + async kernel launch + final sync + D2H.
 * NOT measured: GPU allocation + upload of initial states/constants.
 */
static void bench_stage_aligned_gpu_rk4(benchmark::State& state)
{
    const int num_commuters = static_cast<int>(state.range(0));
    const int G             = static_cast<int>(state.range(1));
    const int NC            = 4 * G;
    const int num_patches   = num_commuters + 1;

    mio::set_log_level(mio::LogLevel::critical);

    if (num_commuters <= 0) {
        state.SkipWithError("num_commuters must be > 0.");
        return;
    }

    // ---- GPU buffers ----
    double *d_mobile = nullptr, *d_lambda_stages = nullptr;
    double *d_rE = nullptr, *d_rI = nullptr;
    int *d_iS = nullptr, *d_iE = nullptr, *d_iI = nullptr, *d_iR = nullptr;

    CUDA_CHECK(cudaMalloc(&d_mobile, (size_t)NC * num_commuters * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lambda_stages, (size_t)4 * G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rE, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rI, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_iS, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iE, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iI, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iR, (size_t)G * sizeof(int)));

    std::vector<double> h_result((size_t)NC * num_commuters);

    for (auto _ : state) {
        for (int patch = 0; patch < num_patches; ++patch) {

            // ====  Setup – NOT measured  ====
            state.PauseTiming();

            ModelType model(G);
            setup_model_gpu(model, static_cast<size_t>(G));

            auto ic = mio::oseir::make_seir_index_cache(model);

            std::vector<double> rate_E(G), rate_I(G);
            for (int g = 0; g < G; ++g) {
                const double tE = model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(g)];
                const double tI = model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(g)];
                rate_E[g]       = (tE > 1e-10) ? (1.0 / tE) : 0.0;
                rate_I[g]       = (tI > 1e-10) ? (1.0 / tI) : 0.0;
            }

            Eigen::VectorXd current_totals = model.populations.get_compartments();
            const double mobile_fraction   = 0.1 * num_commuters;
            const Eigen::VectorXd init_mob = current_totals * (mobile_fraction / num_commuters);

            std::vector<double> h_mobile((size_t)NC * num_commuters);
            for (int cg = 0; cg < num_commuters; ++cg)
                for (int c = 0; c < NC; ++c)
                    h_mobile[(size_t)c * num_commuters + cg] = init_mob[c];

            CUDA_CHECK(cudaMemcpy(d_mobile, h_mobile.data(), (size_t)NC * num_commuters * sizeof(double),
                                  cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_rE, rate_E.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_rI, rate_I.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iS, ic.S.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iE, ic.E.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iI, ic.I.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iR, ic.R.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));

            // Stage caches
            std::vector<double> h_lambda_stages((size_t)4 * G);
            std::vector<Eigen::VectorXd> y_s(4, Eigen::VectorXd(NC));
            std::vector<Eigen::VectorXd> k_s(4, Eigen::VectorXd(NC));
            std::vector<std::vector<double>> lam_s(4, std::vector<double>(G, 0.0));

            auto compute_lambda = [&](const Eigen::VectorXd& y, double t, std::vector<double>& lam) {
                lam.assign(G, 0.0);
                for (int j = 0; j < G; ++j) {
                    const double Nj    = y[ic.S[j]] + y[ic.E[j]] + y[ic.I[j]] + y[ic.R[j]];
                    const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;
                    for (int i = 0; i < G; ++i) {
                        const double coeff =
                            model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>()
                                .get_cont_freq_mat()
                                .get_matrix_at(t)(i, j) *
                            model.parameters
                                .get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(i)] *
                            divNj;
                        lam[i] += coeff * y[ic.I[j]];
                    }
                }
            };

            // ====  Timed loop  ====
            state.ResumeTiming();

            for (ScalarType t = t0_gpu; t < t_max_gpu; t += dt_gpu) {
                if (t + dt_gpu > t_max_gpu + 1e-10)
                    break;

                // CPU: RK4 for aggregate totals
                y_s[0] = current_totals;
                model.get_derivatives(y_s[0], y_s[0], t, k_s[0]);
                compute_lambda(y_s[0], t, lam_s[0]);

                y_s[1] = y_s[0] + (dt_gpu * 0.5) * k_s[0];
                model.get_derivatives(y_s[1], y_s[1], t + 0.5 * dt_gpu, k_s[1]);
                compute_lambda(y_s[1], t + 0.5 * dt_gpu, lam_s[1]);

                y_s[2] = y_s[0] + (dt_gpu * 0.5) * k_s[1];
                model.get_derivatives(y_s[2], y_s[2], t + 0.5 * dt_gpu, k_s[2]);
                compute_lambda(y_s[2], t + 0.5 * dt_gpu, lam_s[2]);

                y_s[3] = y_s[0] + dt_gpu * k_s[2];
                model.get_derivatives(y_s[3], y_s[3], t + dt_gpu, k_s[3]);
                compute_lambda(y_s[3], t + dt_gpu, lam_s[3]);

                Eigen::VectorXd next_totals =
                    current_totals + (dt_gpu / 6.0) * (k_s[0] + 2.0 * k_s[1] + 2.0 * k_s[2] + k_s[3]);

                // Pack 4 × G lambdas for H2D
                for (int s = 0; s < 4; ++s)
                    for (int g = 0; g < G; ++g)
                        h_lambda_stages[(size_t)s * G + g] = lam_s[s][g];

                // H2D: lambda
                CUDA_CHECK(cudaMemcpy(d_lambda_stages, h_lambda_stages.data(), (size_t)4 * G * sizeof(double),
                                      cudaMemcpyHostToDevice));

                // GPU: launch RK4 kernel asynchronously
                launch_seir_commuter_rk4(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, NC, G, dt_gpu,
                                         num_commuters);

                // CPU: advance totals (overlaps with GPU kernel)
                current_totals = next_totals;
            }

            CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(h_result.data(), d_mobile, (size_t)NC * num_commuters * sizeof(double),
                                  cudaMemcpyDeviceToHost));
            benchmark::DoNotOptimize(h_result);
        }
    }

    cudaFree(d_mobile);
    cudaFree(d_lambda_stages);
    cudaFree(d_rE);
    cudaFree(d_rI);
    cudaFree(d_iS);
    cudaFree(d_iE);
    cudaFree(d_iI);
    cudaFree(d_iR);
}

/**
 * @brief GPU stage-aligned Euler benchmark.
 *
 * Simpler than RK4: one derivative call per step, only G lambda values
 * uploaded per step.
 */
static void bench_stage_aligned_gpu_euler(benchmark::State& state)
{
    const int num_commuters = static_cast<int>(state.range(0));
    const int G             = static_cast<int>(state.range(1));
    const int NC            = 4 * G;
    const int num_patches   = num_commuters + 1;

    mio::set_log_level(mio::LogLevel::critical);

    if (num_commuters <= 0) {
        state.SkipWithError("num_commuters must be > 0.");
        return;
    }

    double *d_mobile = nullptr, *d_lambda = nullptr;
    double *d_rE = nullptr, *d_rI = nullptr;
    int *d_iS = nullptr, *d_iE = nullptr, *d_iI = nullptr, *d_iR = nullptr;

    CUDA_CHECK(cudaMalloc(&d_mobile, (size_t)NC * num_commuters * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lambda, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rE, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rI, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_iS, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iE, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iI, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iR, (size_t)G * sizeof(int)));

    std::vector<double> h_result((size_t)NC * num_commuters);

    for (auto _ : state) {
        for (int patch = 0; patch < num_patches; ++patch) {

            state.PauseTiming();

            ModelType model(G);
            setup_model_gpu(model, static_cast<size_t>(G));

            auto ic = mio::oseir::make_seir_index_cache(model);

            std::vector<double> rate_E(G), rate_I(G);
            for (int g = 0; g < G; ++g) {
                const double tE = model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(g)];
                const double tI = model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(g)];
                rate_E[g]       = (tE > 1e-10) ? (1.0 / tE) : 0.0;
                rate_I[g]       = (tI > 1e-10) ? (1.0 / tI) : 0.0;
            }

            Eigen::VectorXd current_totals = model.populations.get_compartments();
            const double mobile_fraction   = 0.1 * num_commuters;
            const Eigen::VectorXd init_mob = current_totals * (mobile_fraction / num_commuters);

            std::vector<double> h_mobile((size_t)NC * num_commuters);
            for (int cg = 0; cg < num_commuters; ++cg)
                for (int c = 0; c < NC; ++c)
                    h_mobile[(size_t)c * num_commuters + cg] = init_mob[c];

            CUDA_CHECK(cudaMemcpy(d_mobile, h_mobile.data(), (size_t)NC * num_commuters * sizeof(double),
                                  cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_rE, rate_E.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_rI, rate_I.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iS, ic.S.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iE, ic.E.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iI, ic.I.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iR, ic.R.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));

            std::vector<double> h_lambda(G);
            Eigen::VectorXd k_tot(NC);

            auto compute_lambda = [&](const Eigen::VectorXd& y, double t, std::vector<double>& lam) {
                lam.assign(G, 0.0);
                for (int j = 0; j < G; ++j) {
                    const double Nj    = y[ic.S[j]] + y[ic.E[j]] + y[ic.I[j]] + y[ic.R[j]];
                    const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;
                    for (int i = 0; i < G; ++i) {
                        const double coeff =
                            model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>()
                                .get_cont_freq_mat()
                                .get_matrix_at(t)(i, j) *
                            model.parameters
                                .get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(i)] *
                            divNj;
                        lam[i] += coeff * y[ic.I[j]];
                    }
                }
            };

            state.ResumeTiming();

            for (ScalarType t = t0_gpu; t < t_max_gpu; t += dt_gpu) {
                if (t + dt_gpu > t_max_gpu + 1e-10)
                    break;

                // CPU: Euler totals + lambda
                model.get_derivatives(current_totals, current_totals, t, k_tot);
                compute_lambda(current_totals, t, h_lambda);
                Eigen::VectorXd next_totals = current_totals + dt_gpu * k_tot;

                // H2D: lambda (G doubles)
                CUDA_CHECK(cudaMemcpy(d_lambda, h_lambda.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));

                // GPU: launch Euler kernel asynchronously
                launch_seir_commuter_euler(d_mobile, d_lambda, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, NC, G, dt_gpu,
                                           num_commuters);

                // CPU: advance totals (overlaps with GPU kernel)
                current_totals = next_totals;
            }

            CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(h_result.data(), d_mobile, (size_t)NC * num_commuters * sizeof(double),
                                  cudaMemcpyDeviceToHost));
            benchmark::DoNotOptimize(h_result);
        }
    }

    cudaFree(d_mobile);
    cudaFree(d_lambda);
    cudaFree(d_rE);
    cudaFree(d_rI);
    cudaFree(d_iS);
    cudaFree(d_iE);
    cudaFree(d_iI);
    cudaFree(d_iR);
}

/**
 * @brief Augmented ODE system [z | Phi_0 | Phi_1 | ... | Phi_{G-1}].
 *
 * Exploits the fact that the fundamental matrix Phi is completely decoupled
 * between age groups -> block-diagonal structure with G independent 4×4 blocks.
 * State vector size: NC + G*16  (vs NC + NC^2 for the full Phi).
 */
struct AugmentedPhiSystemBlockDiag {
    const mio::oseir::Model<ScalarType>& model;
    size_t G;
    size_t NC;

    explicit AugmentedPhiSystemBlockDiag(const mio::oseir::Model<ScalarType>& m)
        : model(m)
        , G(static_cast<size_t>(m.parameters.get_num_groups()))
        , NC(static_cast<size_t>(m.populations.get_num_compartments()))
    {
    }

    void operator()(const Eigen::VectorXd& y, Eigen::VectorXd& dydt, double t)
    {
        const auto z = y.head(static_cast<Eigen::Index>(NC));

        // 1. dz/dt  (standard SEIR RHS)
        Eigen::VectorXd dz(static_cast<Eigen::Index>(NC));
        model.get_derivatives(z, z, t, dz);
        dydt.head(static_cast<Eigen::Index>(NC)) = dz;

        // 2. Force of infection per age group
        using IS            = mio::oseir::InfectionState;
        Eigen::VectorXd foi = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(G));
        for (size_t j = 0; j < G; ++j) {
            const size_t Sj = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Susceptible});
            const size_t Ej = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Exposed});
            const size_t Ij = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Infected});
            const size_t Rj = model.populations.get_flat_index({mio::AgeGroup(static_cast<int>(j)), IS::Recovered});
            const double Nj = z[Sj] + z[Ej] + z[Ij] + z[Rj];
            const double divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;
            for (size_t i = 0; i < G; ++i) {
                const double coeff =
                    model.parameters.template get<mio::oseir::ContactPatterns<ScalarType>>()
                        .get_cont_freq_mat()
                        .get_matrix_at(t)(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) *
                    model.parameters
                        .template get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(
                            static_cast<int>(i))] *
                    divNj;
                foi[static_cast<Eigen::Index>(i)] += coeff * z[Ij];
            }
        }

        // 3. dPhi_g/dt = A_g * Phi_g  for each decoupled 4×4 block
        for (size_t g = 0; g < G; ++g) {
            const double rate_E =
                1.0 / model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];
            const double rate_I =
                1.0 / model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(static_cast<int>(g))];
            const double lambda = foi[static_cast<Eigen::Index>(g)];

            Eigen::Matrix4d Ag = Eigen::Matrix4d::Zero();
            Ag(0, 0)           = -lambda;
            Ag(1, 0)           = lambda;
            Ag(1, 1)           = -rate_E;
            Ag(2, 1)           = rate_E;
            Ag(2, 2)           = -rate_I;
            Ag(3, 2)           = rate_I;

            const auto Phi_g = Eigen::Map<const Eigen::Matrix4d>(y.data() + NC + g * 16);
            Eigen::Map<Eigen::Matrix4d>(dydt.data() + NC + g * 16) = Ag * Phi_g;
        }
    }
};

/**
 * @brief GPU block-diagonal Phi reconstruction benchmark.
 *
 * Algorithm:
 *   1. CPU: integrate augmented ODE (z + G×4×4 Phi blocks) with RK4.
 *      Same work as bench_matrix_phi_reconstruction_blockdiag on the CPU.
 *   2. H2D: upload G×16 Phi doubles
 *   3. GPU: kernel computes Xt[cg] = Phi * X0[cg] for ALL commuter groups
 *      simultaneously (one thread per commuter group, 16 FMAs per age group).
 *   4. cudaDeviceSynchronize + D2H.
 *
 * Measured time: steps 1–4. Steps 1 and 4 are shared with the CPU reference;
 * only step 3 is GPU-specific. The GPU eliminates the O(N×G) matmul loop that
 * dominates the CPU version for large num_commuters.
 */
static void bench_phi_blockdiag_gpu(benchmark::State& state)
{
    const int num_commuters = static_cast<int>(state.range(0));
    const int G             = static_cast<int>(state.range(1));
    const int NC            = 4 * G;
    const int num_patches   = num_commuters + 1;

    mio::set_log_level(mio::LogLevel::critical);

    // ---- Persistent GPU buffers ----
    double *d_X0 = nullptr, *d_Xt = nullptr, *d_Phi = nullptr;
    int *d_iS = nullptr, *d_iE = nullptr, *d_iI = nullptr, *d_iR = nullptr;

    CUDA_CHECK(cudaMalloc(&d_X0, (size_t)NC * num_commuters * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_Xt, (size_t)NC * num_commuters * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_Phi, (size_t)G * 16 * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_iS, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iE, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iI, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iR, (size_t)G * sizeof(int)));

    std::vector<double> h_result((size_t)NC * num_commuters);

    for (auto _ : state) {
        for (int patch = 0; patch < num_patches; ++patch) {

            // ====  Setup – NOT measured  ====
            state.PauseTiming();

            ModelType model(G);
            setup_model_gpu(model, static_cast<size_t>(G));

            auto ic = mio::oseir::make_seir_index_cache(model);

            const Eigen::Index sys_size = NC + G * 16;
            Eigen::VectorXd y(sys_size);
            y.head(NC) = model.populations.get_compartments();
            for (int g = 0; g < G; ++g)
                Eigen::Map<Eigen::Matrix4d>(y.data() + NC + g * 16) = Eigen::Matrix4d::Identity();

            const double mobile_fraction   = 0.1 * num_commuters;
            const Eigen::VectorXd init_mob = model.populations.get_compartments() * (mobile_fraction / num_commuters);

            std::vector<double> h_X0((size_t)NC * num_commuters);
            for (int cg = 0; cg < num_commuters; ++cg)
                for (int c = 0; c < NC; ++c)
                    h_X0[(size_t)c * num_commuters + cg] = init_mob[c];

            CUDA_CHECK(
                cudaMemcpy(d_X0, h_X0.data(), (size_t)NC * num_commuters * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iS, ic.S.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iE, ic.E.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iI, ic.I.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iR, ic.R.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));

            AugmentedPhiSystemBlockDiag sys(model);

            Eigen::VectorXd k1(sys_size), k2(sys_size), k3(sys_size), k4(sys_size);

            // ====  Timed section  ====
            state.ResumeTiming();

            // 1. CPU: integrate augmented ODE over [0, t_max]  (manual RK4)
            double t_cur = t0_gpu;
            while (t_cur < t_max_gpu - 1e-10) {
                const double dt_eff = std::min(static_cast<double>(dt_gpu), t_max_gpu - t_cur);
                sys(y, k1, t_cur);
                sys(y + (dt_eff * 0.5) * k1, k2, t_cur + 0.5 * dt_eff);
                sys(y + (dt_eff * 0.5) * k2, k3, t_cur + 0.5 * dt_eff);
                sys(y + dt_eff * k3, k4, t_cur + dt_eff);
                y += (dt_eff / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                t_cur += dt_eff;
            }

            // 2. Extract Phi matrices and upload to GPU (≤ 2 KB)
            std::vector<double> h_Phi((size_t)G * 16);
            for (int g = 0; g < G; ++g)
                std::copy(y.data() + NC + g * 16, y.data() + NC + g * 16 + 16, h_Phi.data() + g * 16);
            CUDA_CHECK(cudaMemcpy(d_Phi, h_Phi.data(), (size_t)G * 16 * sizeof(double), cudaMemcpyHostToDevice));

            // 3. GPU: Xt[cg] = Phi * X0[cg]  for all commuter groups in parallel
            launch_seir_phi_reconstruction(d_Xt, d_X0, d_Phi, d_iS, d_iE, d_iI, d_iR, NC, G, num_commuters);

            CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(
                cudaMemcpy(h_result.data(), d_Xt, (size_t)NC * num_commuters * sizeof(double), cudaMemcpyDeviceToHost));
            benchmark::DoNotOptimize(h_result);
        }
    }

    cudaFree(d_X0);
    cudaFree(d_Xt);
    cudaFree(d_Phi);
    cudaFree(d_iS);
    cudaFree(d_iE);
    cudaFree(d_iI);
    cudaFree(d_iR);
}

/**
 * @brief Stage-aligned RK4 benchmark with totals ODE also on the GPU.
 *
 * **Patch loop runs on the CPU** – the host iterates sequentially over all
 * P = num_commuters+1 patches.  Per time step, two kernels are launched:
 *   - Kernel 1 <<<1, 1>>>          : totals RK4 for the current patch → d_lambda_stages[4G]
 *   - Kernel 2 <<<ceil(N/256),256>>>: commuter RK4 for N commuters of that patch
 * No H2D transfer occurs inside the timed loop (contrast with bench_stage_aligned_gpu_rk4).
 *
 * @see bench_stage_aligned_allpatches_rk4 – variant where the patch loop is
 *      also moved to the GPU via a 2-D kernel grid.
 */
static void bench_stage_aligned_fullgpu_rk4(benchmark::State& state)
{
    const int num_commuters = static_cast<int>(state.range(0));
    const int G             = static_cast<int>(state.range(1));
    const int NC            = 4 * G;
    const int num_patches   = num_commuters + 1;

    mio::set_log_level(mio::LogLevel::critical);

    if (num_commuters <= 0) {
        state.SkipWithError("num_commuters must be > 0.");
        return;
    }

    // ---- Persistent GPU buffers ----
    double *d_z = nullptr, *d_mobile = nullptr, *d_lambda_stages = nullptr;
    double *d_contact = nullptr, *d_beta = nullptr;
    double *d_rE = nullptr, *d_rI = nullptr;
    int *d_iS = nullptr, *d_iE = nullptr, *d_iI = nullptr, *d_iR = nullptr;

    CUDA_CHECK(cudaMalloc(&d_z, (size_t)NC * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_mobile, (size_t)NC * num_commuters * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lambda_stages, (size_t)4 * G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_contact, (size_t)G * G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_beta, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rE, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rI, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_iS, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iE, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iI, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iR, (size_t)G * sizeof(int)));

    std::vector<double> h_result((size_t)NC * num_commuters);

    for (auto _ : state) {
        for (int patch = 0; patch < num_patches; ++patch) {

            // ====  Setup – NOT measured  ====
            state.PauseTiming();

            ModelType model(G);
            setup_model_gpu(model, static_cast<size_t>(G));

            auto ic = mio::oseir::make_seir_index_cache(model);

            std::vector<double> rate_E(G), rate_I(G), beta(G);
            for (int g = 0; g < G; ++g) {
                const double tE = model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(g)];
                const double tI = model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(g)];
                rate_E[g]       = (tE > 1e-10) ? (1.0 / tE) : 0.0;
                rate_I[g]       = (tI > 1e-10) ? (1.0 / tI) : 0.0;
                beta[g] =
                    model.parameters.get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(g)];
            }

            // Row-major contact matrix [G×G] (contact rate is constant here)
            std::vector<double> h_contact((size_t)G * G);
            for (int i = 0; i < G; ++i)
                for (int j = 0; j < G; ++j)
                    h_contact[(size_t)i * G + j] = model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>()
                                                       .get_cont_freq_mat()
                                                       .get_matrix_at(0.0)(i, j);

            // Initial states
            const Eigen::VectorXd initial_totals = model.populations.get_compartments();
            const double mobile_fraction         = 0.1 * num_commuters;
            const Eigen::VectorXd init_mob       = initial_totals * (mobile_fraction / num_commuters);

            std::vector<double> h_z(NC);
            for (int c = 0; c < NC; ++c)
                h_z[c] = initial_totals[c];

            std::vector<double> h_mobile((size_t)NC * num_commuters);
            for (int cg = 0; cg < num_commuters; ++cg)
                for (int c = 0; c < NC; ++c)
                    h_mobile[(size_t)c * num_commuters + cg] = init_mob[c];

            CUDA_CHECK(cudaMemcpy(d_z, h_z.data(), (size_t)NC * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_mobile, h_mobile.data(), (size_t)NC * num_commuters * sizeof(double),
                                  cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_contact, h_contact.data(), (size_t)G * G * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_beta, beta.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_rE, rate_E.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_rI, rate_I.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iS, ic.S.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iE, ic.E.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iI, ic.I.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(d_iR, ic.R.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));

            // ====  Timed section  ====
            state.ResumeTiming();

            for (ScalarType t = t0_gpu; t < t_max_gpu; t += dt_gpu) {
                if (t + dt_gpu > t_max_gpu + 1e-10)
                    break;

                // Kernel 1 (1 thread): totals RK4 + write d_lambda_stages
                launch_seir_totals_rk4_lambda(d_z, d_lambda_stages, d_contact, d_beta, d_rE, d_rI, d_iS, d_iE, d_iI,
                                              d_iR, NC, G, dt_gpu);

                // Kernel 2 (N threads): commuter RK4
                launch_seir_commuter_rk4(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, NC, G, dt_gpu,
                                         num_commuters);
            }

            CUDA_CHECK(cudaDeviceSynchronize());
            CUDA_CHECK(cudaMemcpy(h_result.data(), d_mobile, (size_t)NC * num_commuters * sizeof(double),
                                  cudaMemcpyDeviceToHost));
            benchmark::DoNotOptimize(h_result);
        }
    }

    cudaFree(d_z);
    cudaFree(d_mobile);
    cudaFree(d_lambda_stages);
    cudaFree(d_contact);
    cudaFree(d_beta);
    cudaFree(d_rE);
    cudaFree(d_rI);
    cudaFree(d_iS);
    cudaFree(d_iE);
    cudaFree(d_iI);
    cudaFree(d_iR);
}

/**
 * @brief Stage-aligned RK4 with the outer patch loop removed from the CPU.
 *
 * **Patch loop runs on the GPU** – all P = num_commuters+1 patches are
 * integrated simultaneously in every time step.  Per time step, two kernels
 * are launched once for *all* patches:
 *   - Kernel 1 <<<ceil(P/256), 256>>>         : totals RK4 for all P patches → d_lambda_stages[4G×P]
 *   - Kernel 2 <<<dim3(ceil(N/256),P), 256>>>: commuter RK4 for all P×N pairs
 * No H2D transfer and no CPU patch loop in the timed section.
 * Trade-off: maximises GPU parallelism but requires O(P·N) device memory.
 *
 * @see bench_stage_aligned_fullgpu_rk4 – variant where the patch loop is
 *      still on the CPU (sequential over patches, parallel over commuters).
 */
static void bench_stage_aligned_allpatches_rk4(benchmark::State& state)
{
    const int N  = static_cast<int>(state.range(0)); // commuters per patch
    const int G  = static_cast<int>(state.range(1)); // age groups
    const int P  = N + 1; // number of patches
    const int NC = 4 * G;

    mio::set_log_level(mio::LogLevel::critical);

    if (N <= 0 || P <= 0) {
        state.SkipWithError("N and P must be > 0.");
        return;
    }

    // ---- GPU buffers (allocated once before the benchmark loop) ----
    double *d_z = nullptr, *d_mobile = nullptr, *d_lambda_stages = nullptr;
    double *d_contact = nullptr, *d_beta = nullptr;
    double *d_rE = nullptr, *d_rI = nullptr;
    int *d_iS = nullptr, *d_iE = nullptr, *d_iI = nullptr, *d_iR = nullptr;

    CUDA_CHECK(cudaMalloc(&d_z, (size_t)NC * P * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_lambda_stages, (size_t)4 * G * P * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_mobile, (size_t)NC * P * N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_contact, (size_t)G * G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_beta, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rE, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_rI, (size_t)G * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_iS, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iE, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iI, (size_t)G * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_iR, (size_t)G * sizeof(int)));

    std::vector<double> h_z((size_t)NC * P);
    std::vector<double> h_mobile((size_t)NC * P * N);

    for (auto _ : state) {

        // ====  Setup – NOT measured  ====
        state.PauseTiming();

        ModelType model(G);
        setup_model_gpu(model, static_cast<size_t>(G));

        auto ic = mio::oseir::make_seir_index_cache(model);

        std::vector<double> rate_E(G), rate_I(G), beta(G);
        for (int g = 0; g < G; ++g) {
            const double tE = model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(g)];
            const double tI = model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(g)];
            rate_E[g]       = (tE > 1e-10) ? (1.0 / tE) : 0.0;
            rate_I[g]       = (tI > 1e-10) ? (1.0 / tI) : 0.0;
            beta[g] =
                model.parameters.get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(g)];
        }

        std::vector<double> h_contact((size_t)G * G);
        for (int i = 0; i < G; ++i)
            for (int j = 0; j < G; ++j)
                h_contact[(size_t)i * G + j] =
                    model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>().get_cont_freq_mat().get_matrix_at(
                        0.0)(i, j);

        // Initial condition (same for every patch)
        const Eigen::VectorXd initial_totals = model.populations.get_compartments();
        const double mob_frac                = 0.1 * N;
        const Eigen::VectorXd init_mob       = initial_totals * (mob_frac / N);

        for (int c = 0; c < NC; ++c)
            for (int p = 0; p < P; ++p)
                h_z[(size_t)c * P + p] = initial_totals[c];

        for (int p = 0; p < P; ++p)
            for (int c = 0; c < NC; ++c)
                std::fill(h_mobile.begin() + (size_t)(p * NC + c) * N, h_mobile.begin() + (size_t)(p * NC + c) * N + N,
                          init_mob[c]);

        CUDA_CHECK(cudaMemcpy(d_z, h_z.data(), (size_t)NC * P * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_mobile, h_mobile.data(), (size_t)NC * P * N * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_contact, h_contact.data(), (size_t)G * G * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_beta, beta.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_rE, rate_E.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_rI, rate_I.data(), (size_t)G * sizeof(double), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_iS, ic.S.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_iE, ic.E.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_iI, ic.I.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_iR, ic.R.data(), (size_t)G * sizeof(int), cudaMemcpyHostToDevice));

        // ====  Timed section  ====
        state.ResumeTiming();

        for (ScalarType t = t0_gpu; t < t_max_gpu; t += dt_gpu) {
            if (t + dt_gpu > t_max_gpu + 1e-10)
                break;

            // Kernel 1: totals RK4 for all P patches (<<<ceil(P/256), 256>>>)
            launch_seir_totals_rk4_lambda_allpatches(d_z, d_lambda_stages, d_contact, d_beta, d_rE, d_rI, d_iS, d_iE,
                                                     d_iI, d_iR, NC, G, P, dt_gpu);

            // Kernel 2: commuter RK4 for all P×N pairs (<<<dim3(ceil(N/256),P), 256>>>)
            // Automatically ordered after kernel 1 (same default stream).
            launch_seir_commuter_rk4_allpatches(d_mobile, d_lambda_stages, d_rE, d_rI, d_iS, d_iE, d_iI, d_iR, NC, G, P,
                                                N, dt_gpu);
        }

        // Barrier – ensures all GPU work is done before the timer stops.
        CUDA_CHECK(cudaDeviceSynchronize());
        benchmark::DoNotOptimize(d_mobile);
    }

    cudaFree(d_z);
    cudaFree(d_lambda_stages);
    cudaFree(d_mobile);
    cudaFree(d_contact);
    cudaFree(d_beta);
    cudaFree(d_rE);
    cudaFree(d_rI);
    cudaFree(d_iS);
    cudaFree(d_iE);
    cudaFree(d_iI);
    cudaFree(d_iR);
}

} // namespace benchmark_mio
} // namespace mio

// ---------------------------------------------------------------------------
// Benchmark registration
// ---------------------------------------------------------------------------

BENCHMARK(mio::benchmark_mio::bench_stage_aligned_gpu_rk4)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts_gpu)
            for (int g : mio::benchmark_mio::age_group_counts_gpu)
                b->Args({i, g});
    }) -> Name("stage-aligned-gpu(RK4)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_stage_aligned_gpu_euler)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts_gpu)
            for (int g : mio::benchmark_mio::age_group_counts_gpu)
                b->Args({i, g});
    }) -> Name("stage-aligned-gpu(Euler)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_phi_blockdiag_gpu)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts_gpu)
            for (int g : mio::benchmark_mio::age_group_counts_gpu)
                b->Args({i, g});
    }) -> Name("phi-blockdiag-gpu") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_stage_aligned_fullgpu_rk4)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts_gpu)
            for (int g : mio::benchmark_mio::age_group_counts_gpu)
                b->Args({i, g});
    }) -> Name("stage-aligned-fullgpu(RK4)") -> Unit(::benchmark::kMicrosecond);

// All-patches: outer patch loop runs on the GPU via a 2D kernel grid.
// P = num_commuters + 1 (same as bench_stage_aligned_fullgpu(RK4) for direct comparison).
BENCHMARK(mio::benchmark_mio::bench_stage_aligned_allpatches_rk4)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts_gpu)
            for (int g : mio::benchmark_mio::age_group_counts_gpu)
                b->Args({i, g});
    }) -> Name("stage-aligned-allpatches(RK4)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK_MAIN();
