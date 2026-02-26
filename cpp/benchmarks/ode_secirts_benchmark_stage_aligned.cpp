#include "benchmark/benchmark.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/math/euler.h"

#include "ode_secirts/model.h"
#include "ode_secirts/parameters.h"
#include "ode_secirts/state_estimators.h"

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <vector>
#include <iostream>

using ScalarType = double;

namespace mio
{
namespace benchmark_mio
{

const ScalarType t0    = 0.0;
const ScalarType t_max = 0.5;
const ScalarType dt    = 0.1;

const ScalarType t_max_phi = 0.5;
const ScalarType dt_phi    = 0.1;

namespace
{
struct BenchSetupPrinterSecirts {
    BenchSetupPrinterSecirts()
    {
        std::cout << "SECIRTS benchmark setup:\n";
        std::cout << "  t_max = " << t_max << "\n";
        std::cout << "  dt    = " << dt << "\n";
        std::cout << "  t_max_phi = " << t_max_phi << "\n";
        std::cout << "  dt_phi    = " << dt_phi << "\n";
        std::cout << std::flush;
    }
};
static BenchSetupPrinterSecirts bench_setup_printer_secirts;
} // namespace

const std::vector<int> commuter_group_counts = {16, 32, 64, 128, 256, 512, 1024};
const std::vector<int> age_group_counts      = {1, 2, 3, 4, 5, 6};

template <class ModelType>
void setup_parameters(ModelType& model, size_t num_age_groups)
{
    for (size_t g = 0; g < num_age_groups; ++g) {
        auto ag                                                                                    = mio::AgeGroup(g);
        model.parameters.template get<mio::osecirts::TimeExposed<ScalarType>>()[ag]                = 3.33;
        model.parameters.template get<mio::osecirts::TimeInfectedNoSymptoms<ScalarType>>()[ag]     = 1.87;
        model.parameters.template get<mio::osecirts::TimeInfectedSymptoms<ScalarType>>()[ag]       = 7;
        model.parameters.template get<mio::osecirts::TimeInfectedSevere<ScalarType>>()[ag]         = 6;
        model.parameters.template get<mio::osecirts::TimeInfectedCritical<ScalarType>>()[ag]       = 7;
        model.parameters.template get<mio::osecirts::TimeTemporaryImmunityPI<ScalarType>>()[ag]    = 60;
        model.parameters.template get<mio::osecirts::TimeTemporaryImmunityII<ScalarType>>()[ag]    = 60;
        model.parameters.template get<mio::osecirts::TimeWaningPartialImmunity<ScalarType>>()[ag]  = 180;
        model.parameters.template get<mio::osecirts::TimeWaningImprovedImmunity<ScalarType>>()[ag] = 180;

        model.parameters.template get<mio::osecirts::TransmissionProbabilityOnContact<ScalarType>>()[ag]  = 0.15;
        model.parameters.template get<mio::osecirts::RelativeTransmissionNoSymptoms<ScalarType>>()[ag]    = 0.5;
        model.parameters.template get<mio::osecirts::RiskOfInfectionFromSymptomatic<ScalarType>>()[ag]    = 0.0;
        model.parameters.template get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<ScalarType>>()[ag] = 0.4;
        model.parameters.template get<mio::osecirts::RecoveredPerInfectedNoSymptoms<ScalarType>>()[ag]    = 0.2;
        model.parameters.template get<mio::osecirts::SeverePerInfectedSymptoms<ScalarType>>()[ag]         = 0.1;
        model.parameters.template get<mio::osecirts::CriticalPerSevere<ScalarType>>()[ag]                 = 0.1;
        model.parameters.template get<mio::osecirts::DeathsPerCritical<ScalarType>>()[ag]                 = 0.1;

        model.parameters.template get<mio::osecirts::ReducExposedPartialImmunity<ScalarType>>()[ag]           = 0.8;
        model.parameters.template get<mio::osecirts::ReducExposedImprovedImmunity<ScalarType>>()[ag]          = 0.331;
        model.parameters.template get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<ScalarType>>()[ag]  = 0.65;
        model.parameters.template get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[ag] = 0.243;
        model.parameters.template get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[ag] =
            0.1;
        model.parameters
            .template get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[ag] = 0.091;
        model.parameters.template get<mio::osecirts::ReducTimeInfectedMild<ScalarType>>()[ag]               = 0.9;
    }

    model.parameters.template get<mio::osecirts::ICUCapacity<ScalarType>>()          = 100;
    model.parameters.template get<mio::osecirts::TestAndTraceCapacity<ScalarType>>() = 0.0143;

    model.parameters.template get<mio::osecirts::DailyPartialVaccinations<ScalarType>>().resize(mio::SimulationDay(10));
    model.parameters.template get<mio::osecirts::DailyFullVaccinations<ScalarType>>().resize(mio::SimulationDay(10));
    model.parameters.template get<mio::osecirts::DailyBoosterVaccinations<ScalarType>>().resize(mio::SimulationDay(10));

    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.template get<mio::osecirts::ContactPatterns<ScalarType>>();
    const double fact = 1.0 / (double)num_age_groups;
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_age_groups, num_age_groups, fact * 10.0));
    model.parameters.template get<mio::osecirts::Seasonality<ScalarType>>() = 0.2;
}

std::vector<double> get_base_population()
{
    using IS = mio::osecirts::InfectionState;
    std::vector<double> base_pop((size_t)IS::Count, 0.0);
    base_pop[(size_t)IS::ExposedNaive]                       = 20;
    base_pop[(size_t)IS::ExposedImprovedImmunity]            = 20;
    base_pop[(size_t)IS::ExposedPartialImmunity]             = 20;
    base_pop[(size_t)IS::InfectedNoSymptomsNaive]            = 30;
    base_pop[(size_t)IS::InfectedNoSymptomsPartialImmunity]  = 30;
    base_pop[(size_t)IS::InfectedNoSymptomsImprovedImmunity] = 30;
    base_pop[(size_t)IS::InfectedSymptomsNaive]              = 40;
    base_pop[(size_t)IS::InfectedSymptomsPartialImmunity]    = 40;
    base_pop[(size_t)IS::InfectedSymptomsImprovedImmunity]   = 40;
    base_pop[(size_t)IS::InfectedSevereNaive]                = 30;
    base_pop[(size_t)IS::InfectedSevereImprovedImmunity]     = 30;
    base_pop[(size_t)IS::InfectedSeverePartialImmunity]      = 30;
    base_pop[(size_t)IS::InfectedCriticalNaive]              = 20;
    base_pop[(size_t)IS::InfectedCriticalPartialImmunity]    = 20;
    base_pop[(size_t)IS::InfectedCriticalImprovedImmunity]   = 20;
    base_pop[(size_t)IS::SusceptibleNaive]                   = 1000;
    base_pop[(size_t)IS::SusceptiblePartialImmunity]         = 1200;
    base_pop[(size_t)IS::SusceptibleImprovedImmunity]        = 1000;
    base_pop[(size_t)IS::TemporaryImmunePartialImmunity]     = 60;
    base_pop[(size_t)IS::TemporaryImmuneImprovedImmunity]    = 70;
    return base_pop;
}

// ==============================================================================
// 1. STAGE ALIGNED RK4
// ==============================================================================
static void bench_stage_aligned_rk4(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    struct SecirtsRK4StageCache {
        std::vector<Eigen::VectorXd> y;
        std::vector<Eigen::VectorXd> k_tot;
        std::vector<mio::osecirts::SecirtsIntensities> intensities;

        void resize(size_t nc, size_t num_age_groups)
        {
            y.resize(4, Eigen::VectorXd(nc));
            k_tot.resize(4, Eigen::VectorXd(nc));
            intensities.resize(4);
            for (int i = 0; i < 4; ++i) {
                intensities[i].resize(num_age_groups);
            }
        }
    };

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();

            mio::osecirts::Model<ScalarType> model(num_age_groups);
            setup_parameters(model, num_age_groups);

            auto base_pop        = get_base_population();
            double frac_commuter = 0.3 / num_commuter_groups;

            Eigen::VectorXd current_totals(num_age_groups * 29);
            std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, Eigen::VectorXd(num_age_groups * 29));

            for (size_t g = 0; g < num_age_groups; ++g) {
                for (size_t state_idx = 0; state_idx < 29; ++state_idx) {
                    model.populations[{mio::AgeGroup(g), (mio::osecirts::InfectionState)state_idx}] =
                        base_pop[state_idx];
                    size_t flat_idx =
                        model.populations.get_flat_index({mio::AgeGroup(g), (mio::osecirts::InfectionState)state_idx});
                    current_totals[flat_idx] = base_pop[state_idx];
                    for (int cg = 0; cg < num_commuter_groups; ++cg) {
                        mobile_pops[cg][flat_idx] = base_pop[state_idx] * frac_commuter;
                    }
                }
            }

            mio::osecirts::ModelEvaluatorSecirts eval(model);
            size_t NC = eval.NC;
            SecirtsRK4StageCache cache;
            cache.resize(NC, num_age_groups);

            Eigen::VectorXd k1_com(NC), k2_com(NC), k3_com(NC), k4_com(NC);
            Eigen::VectorXd c2(NC), c3(NC), c4(NC);

            state.ResumeTiming();

            for (ScalarType t = t0; t < t_max; t += dt) {
                if (t + dt > t_max + 1e-10)
                    break;

                // calculate RK4 stages for totals
                // Stage 1
                cache.y[0] = current_totals;
                eval.get_totals_rhs(cache.y[0], t, cache.k_tot[0]);
                eval.get_intensities(cache.y[0], t, cache.intensities[0]);

                // Stage 2
                cache.y[1] = cache.y[0] + (dt * 0.5) * cache.k_tot[0];
                eval.get_totals_rhs(cache.y[1], t + 0.5 * dt, cache.k_tot[1]);
                eval.get_intensities(cache.y[1], t + 0.5 * dt, cache.intensities[1]);

                // Stage 3
                cache.y[2] = cache.y[0] + (dt * 0.5) * cache.k_tot[1];
                eval.get_totals_rhs(cache.y[2], t + 0.5 * dt, cache.k_tot[2]);
                eval.get_intensities(cache.y[2], t + 0.5 * dt, cache.intensities[2]);

                // Stage 4
                cache.y[3] = cache.y[0] + dt * cache.k_tot[2];
                eval.get_totals_rhs(cache.y[3], t + dt, cache.k_tot[3]);
                eval.get_intensities(cache.y[3], t + dt, cache.intensities[3]);

                Eigen::VectorXd next_totals = current_totals + (dt / 6.0) * (cache.k_tot[0] + 2.0 * cache.k_tot[1] +
                                                                             2.0 * cache.k_tot[2] + cache.k_tot[3]);

                // reconstruct
                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    Eigen::Ref<Eigen::VectorXd> Xc = mobile_pops[cg];

                    // Stage 1
                    eval.get_commuter_rhs(cache.intensities[0], Xc, k1_com);

                    // Stage 2
                    c2 = Xc + (dt * 0.5) * k1_com;
                    eval.get_commuter_rhs(cache.intensities[1], c2, k2_com);

                    // Stage 3
                    c3 = Xc + (dt * 0.5) * k2_com;
                    eval.get_commuter_rhs(cache.intensities[2], c3, k3_com);

                    // Stage 4
                    c4 = Xc + dt * k3_com;
                    eval.get_commuter_rhs(cache.intensities[3], c4, k4_com);

                    // Final update
                    Xc += (dt / 6.0) * (k1_com + 2.0 * k2_com + 2.0 * k3_com + k4_com);
                }

                current_totals = next_totals;
            }
            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

// ==============================================================================
// 2. STAGE ALIGNED EULER
// ==============================================================================
static void bench_stage_aligned_euler(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;

    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();

            mio::osecirts::Model<ScalarType> model(num_age_groups);
            setup_parameters(model, num_age_groups);

            auto base_pop        = get_base_population();
            double frac_commuter = 0.3 / num_commuter_groups;

            Eigen::VectorXd current_totals(num_age_groups * 29);
            std::vector<Eigen::VectorXd> mobile_pops(num_commuter_groups, Eigen::VectorXd(num_age_groups * 29));

            for (size_t g = 0; g < num_age_groups; ++g) {
                for (size_t state_idx = 0; state_idx < 29; ++state_idx) {
                    model.populations[{mio::AgeGroup(g), (mio::osecirts::InfectionState)state_idx}] =
                        base_pop[state_idx];
                    size_t flat_idx =
                        model.populations.get_flat_index({mio::AgeGroup(g), (mio::osecirts::InfectionState)state_idx});
                    current_totals[flat_idx] = base_pop[state_idx];
                    for (int cg = 0; cg < num_commuter_groups; ++cg) {
                        mobile_pops[cg][flat_idx] = base_pop[state_idx] * frac_commuter;
                    }
                }
            }

            mio::osecirts::ModelEvaluatorSecirts eval(model);
            size_t NC = eval.NC;

            Eigen::VectorXd k1_tot(NC);
            Eigen::VectorXd k1_com(NC);
            mio::osecirts::SecirtsIntensities int1;

            state.ResumeTiming();

            for (ScalarType t = t0; t < t_max; t += dt) {
                if (t + dt > t_max + 1e-10)
                    break;

                eval.get_totals_rhs(current_totals, t, k1_tot);
                eval.get_intensities(current_totals, t, int1);

                for (int cg = 0; cg < num_commuter_groups; ++cg) {
                    eval.get_commuter_rhs(int1, mobile_pops[cg], k1_com);
                    mobile_pops[cg] += dt * k1_com;
                }
                current_totals += dt * k1_tot;
            }
            benchmark::DoNotOptimize(mobile_pops);
        }
    }
}

// ==============================================================================
// 3. STANDARD LAGRANGIAN RK4
// ==============================================================================
static void bench_standard_lagrangian_rk4(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();

            mio::osecirts::ModelSecirtsExplicit<ScalarType> model_explicit(num_age_groups, num_commuter_groups);
            setup_parameters(model_explicit, num_age_groups);

            auto base_pop            = get_base_population();
            double frac_non_commuter = 0.7;
            double frac_commuter     = 0.3 / num_commuter_groups;

            for (size_t g = 0; g < num_age_groups; ++g) {
                for (size_t state_idx = 0; state_idx < 29; ++state_idx) {
                    model_explicit.populations[{mio::AgeGroup(g), mio::osecirts::CommuterType::NonCommuter,
                                                (mio::osecirts::InfectionState)state_idx}] =
                        base_pop[state_idx] * frac_non_commuter;
                    for (int c = 0; c < num_commuter_groups; ++c) {
                        auto c_type = static_cast<mio::osecirts::CommuterType>(
                            static_cast<int>(mio::osecirts::CommuterType::CommuterBase) + c);
                        model_explicit
                            .populations[{mio::AgeGroup(g), c_type, (mio::osecirts::InfectionState)state_idx}] =
                            base_pop[state_idx] * frac_commuter;
                    }
                }
            }

            mio::FlowSimulation<ScalarType, mio::osecirts::ModelSecirtsExplicit<ScalarType>> sim_explicit(
                model_explicit, t0, dt);
            auto integrator_rk =
                std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
            sim_explicit.set_integrator(integrator_rk);

            state.ResumeTiming();
            sim_explicit.advance(t_max);
        }
    }
}

// ==============================================================================
// 4. STANDARD LAGRANGIAN EULER
// ==============================================================================
static void bench_standard_lagrangian_euler(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {

            state.PauseTiming();

            mio::osecirts::ModelSecirtsExplicit<ScalarType> model_explicit(num_age_groups, num_commuter_groups);
            setup_parameters(model_explicit, num_age_groups);

            auto base_pop            = get_base_population();
            double frac_non_commuter = 0.7;
            double frac_commuter     = 0.3 / num_commuter_groups;

            for (size_t g = 0; g < num_age_groups; ++g) {
                for (size_t state_idx = 0; state_idx < 29; ++state_idx) {
                    model_explicit.populations[{mio::AgeGroup(g), mio::osecirts::CommuterType::NonCommuter,
                                                (mio::osecirts::InfectionState)state_idx}] =
                        base_pop[state_idx] * frac_non_commuter;
                    for (int c = 0; c < num_commuter_groups; ++c) {
                        auto c_type = static_cast<mio::osecirts::CommuterType>(
                            static_cast<int>(mio::osecirts::CommuterType::CommuterBase) + c);
                        model_explicit
                            .populations[{mio::AgeGroup(g), c_type, (mio::osecirts::InfectionState)state_idx}] =
                            base_pop[state_idx] * frac_commuter;
                    }
                }
            }

            mio::FlowSimulation<ScalarType, mio::osecirts::ModelSecirtsExplicit<ScalarType>> sim_explicit(
                model_explicit, t0, dt);
            auto integrator_euler = std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
            sim_explicit.set_integrator(integrator_euler);

            state.ResumeTiming();
            sim_explicit.advance(t_max);
        }
    }
}

// ==============================================================================
// 5. MATRIX PHI RECONSTRUCTION
// ==============================================================================
template <class Integrator>
static void bench_matrix_phi_reconstruction(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();

            mio::osecirts::Model<ScalarType> model(num_age_groups);
            setup_parameters(model, num_age_groups);

            auto base_pop        = get_base_population();
            double frac_commuter = 0.3 / num_commuter_groups;

            size_t NC       = 29 * num_age_groups;
            size_t sys_size = NC + NC * NC;

            Eigen::VectorXd initial_totals = Eigen::VectorXd::Zero(NC);
            Eigen::MatrixXd X0             = Eigen::MatrixXd::Zero(NC, num_commuter_groups);

            for (size_t g = 0; g < num_age_groups; ++g) {
                for (size_t state_idx = 0; state_idx < 29; ++state_idx) {
                    model.populations[{mio::AgeGroup(g), (mio::osecirts::InfectionState)state_idx}] =
                        base_pop[state_idx];
                    size_t flat_idx =
                        model.populations.get_flat_index({mio::AgeGroup(g), (mio::osecirts::InfectionState)state_idx});
                    initial_totals[flat_idx] = base_pop[state_idx];
                    for (int cg = 0; cg < num_commuter_groups; ++cg) {
                        X0(flat_idx, cg) = base_pop[state_idx] * frac_commuter;
                    }
                }
            }

            Eigen::VectorXd y0(sys_size);
            y0.head(NC)        = initial_totals;
            Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(NC, NC);
            y0.tail(NC * NC)   = Eigen::Map<Eigen::VectorXd>(Id.data(), NC * NC);

            Integrator stepper;
            mio::osecirts::AugmentedPhiSystemSecirts sys(model);

            Eigen::MatrixXd Xt(NC, num_commuter_groups);

            state.ResumeTiming();

            double t          = t0;
            Eigen::VectorXd y = y0;
            while (t < t_max_phi - 1e-10) {
                double dt_eff = std::min(dt_phi, t_max_phi - t);
                stepper.do_step(sys, y, t, dt_eff);
                t += dt_eff;
            }

            const auto Phi_final = Eigen::Map<const Eigen::MatrixXd>(y.tail(NC * NC).data(), NC, NC);
            Xt.noalias()         = Phi_final * X0;

            benchmark::DoNotOptimize(Xt);
        }
    }
}

// ==============================================================================
// 6. MATRIX PHI RECONSTRUCTION (Block-Diagonal)
// ==============================================================================
template <class Integrator>
static void bench_matrix_phi_reconstruction_blockdiag(::benchmark::State& state)
{
    const int num_commuter_groups = state.range(0);
    const size_t num_age_groups   = static_cast<size_t>(state.range(1));
    const auto num_patches        = num_commuter_groups + 1;
    mio::set_log_level(mio::LogLevel::critical);

    for (auto _ : state) {
        for (auto patch = 0; patch < num_patches; patch++) {
            state.PauseTiming();

            mio::osecirts::Model<ScalarType> model(num_age_groups);
            setup_parameters(model, num_age_groups);

            auto base_pop        = get_base_population();
            double frac_commuter = 0.3 / num_commuter_groups;

            size_t G              = num_age_groups;
            size_t block_size     = 29;
            size_t block_elements = block_size * block_size;
            size_t NC             = block_size * G;
            size_t sys_size       = NC + G * block_elements;

            Eigen::VectorXd initial_totals = Eigen::VectorXd::Zero(NC);
            Eigen::MatrixXd X0             = Eigen::MatrixXd::Zero(NC, num_commuter_groups);

            for (size_t g = 0; g < G; ++g) {
                for (size_t state_idx = 0; state_idx < block_size; ++state_idx) {
                    model.populations[{mio::AgeGroup(g), (mio::osecirts::InfectionState)state_idx}] =
                        base_pop[state_idx];
                    size_t flat_idx =
                        model.populations.get_flat_index({mio::AgeGroup(g), (mio::osecirts::InfectionState)state_idx});
                    initial_totals[flat_idx] = base_pop[state_idx];
                    for (int cg = 0; cg < num_commuter_groups; ++cg) {
                        X0(flat_idx, cg) = base_pop[state_idx] * frac_commuter;
                    }
                }
            }

            Eigen::VectorXd y0(sys_size);
            y0.head(NC) = initial_totals;
            for (size_t g = 0; g < G; ++g) {
                Eigen::Map<Eigen::MatrixXd>(y0.data() + NC + g * block_elements, block_size, block_size) =
                    Eigen::MatrixXd::Identity(block_size, block_size);
            }

            Integrator stepper;
            mio::osecirts::AugmentedPhiSystemSecirtsBlockDiag sys(model);

            Eigen::MatrixXd Xt(NC, num_commuter_groups);

            state.ResumeTiming();

            double t          = t0;
            Eigen::VectorXd y = y0;
            while (t < t_max_phi - 1e-10) {
                double dt_eff = std::min(dt_phi, t_max_phi - t);
                stepper.do_step(sys, y, t, dt_eff);
                t += dt_eff;
            }

            for (size_t g = 0; g < G; ++g) {
                const auto Phi_g =
                    Eigen::Map<const Eigen::MatrixXd>(y.data() + NC + g * block_elements, block_size, block_size);
                Xt.block(block_size * g, 0, block_size, num_commuter_groups).noalias() =
                    Phi_g * X0.block(block_size * g, 0, block_size, num_commuter_groups);
            }

            benchmark::DoNotOptimize(Xt);
        }
    }
}

static void bench_matrix_phi_rk4(benchmark::State& state)
{
    bench_matrix_phi_reconstruction<boost::numeric::odeint::runge_kutta4<
        Eigen::VectorXd, double, Eigen::VectorXd, double, boost::numeric::odeint::vector_space_algebra>>(state);
}

static void bench_matrix_phi_euler(benchmark::State& state)
{
    bench_matrix_phi_reconstruction<boost::numeric::odeint::euler<Eigen::VectorXd, double, Eigen::VectorXd, double,
                                                                  boost::numeric::odeint::vector_space_algebra>>(state);
}

static void bench_matrix_phi_blockdiag_rk4(benchmark::State& state)
{
    bench_matrix_phi_reconstruction_blockdiag<boost::numeric::odeint::runge_kutta4<
        Eigen::VectorXd, double, Eigen::VectorXd, double, boost::numeric::odeint::vector_space_algebra>>(state);
}

static void bench_matrix_phi_blockdiag_euler(benchmark::State& state)
{
    bench_matrix_phi_reconstruction_blockdiag<boost::numeric::odeint::euler<
        Eigen::VectorXd, double, Eigen::VectorXd, double, boost::numeric::odeint::vector_space_algebra>>(state);
}

} // namespace benchmark_mio
} // namespace mio

BENCHMARK(mio::benchmark_mio::bench_stage_aligned_rk4)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts) {
            for (int g : mio::benchmark_mio::age_group_counts) {
                b->Args({i, g});
            }
        }
    }) -> Name("stage-aligned(RK4)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_stage_aligned_euler)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts) {
            for (int g : mio::benchmark_mio::age_group_counts) {
                b->Args({i, g});
            }
        }
    }) -> Name("stage-aligned(Euler)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_standard_lagrangian_rk4)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts) {
            for (int g : mio::benchmark_mio::age_group_counts) {
                b->Args({i, g});
            }
        }
    }) -> Name("lagrange_rk4") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_standard_lagrangian_euler)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts) {
            for (int g : mio::benchmark_mio::age_group_counts) {
                b->Args({i, g});
            }
        }
    }) -> Name("lagrange_euler") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_matrix_phi_rk4)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts) {
            for (int g : mio::benchmark_mio::age_group_counts) {
                b->Args({i, g});
            }
        }
    }) -> Name("matrix_phi(RK4)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_matrix_phi_euler)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts) {
            for (int g : mio::benchmark_mio::age_group_counts) {
                b->Args({i, g});
            }
        }
    }) -> Name("matrix_phi(Euler)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_matrix_phi_blockdiag_rk4)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts) {
            for (int g : mio::benchmark_mio::age_group_counts) {
                b->Args({i, g});
            }
        }
    }) -> Name("matrix_phi_blockdiag(RK4)") -> Unit(::benchmark::kMicrosecond);

BENCHMARK(mio::benchmark_mio::bench_matrix_phi_blockdiag_euler)
    ->Apply([](auto* b) {
        for (int i : mio::benchmark_mio::commuter_group_counts) {
            for (int g : mio::benchmark_mio::age_group_counts) {
                b->Args({i, g});
            }
        }
    }) -> Name("matrix_phi_blockdiag(Euler)") -> Unit(::benchmark::kMicrosecond);

// run all benchmarks
BENCHMARK_MAIN();