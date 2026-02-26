#include "ode_secirts/model.h"
#include "ode_secirts/parameters.h"
#include "ode_secirts/state_estimators.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/math/stepper_wrapper.h"
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>

using namespace mio::osecirts;

template <class ModelType>
void setup_parameters(ModelType& model)
{
    model.parameters.template get<TimeExposed<double>>()[mio::AgeGroup(0)]                = 3.33;
    model.parameters.template get<TimeInfectedNoSymptoms<double>>()[mio::AgeGroup(0)]     = 1.87;
    model.parameters.template get<TimeInfectedSymptoms<double>>()[mio::AgeGroup(0)]       = 7;
    model.parameters.template get<TimeInfectedSevere<double>>()[mio::AgeGroup(0)]         = 6;
    model.parameters.template get<TimeInfectedCritical<double>>()[mio::AgeGroup(0)]       = 7;
    model.parameters.template get<TimeTemporaryImmunityPI<double>>()[mio::AgeGroup(0)]    = 60;
    model.parameters.template get<TimeTemporaryImmunityII<double>>()[mio::AgeGroup(0)]    = 60;
    model.parameters.template get<TimeWaningPartialImmunity<double>>()[mio::AgeGroup(0)]  = 180;
    model.parameters.template get<TimeWaningImprovedImmunity<double>>()[mio::AgeGroup(0)] = 180;

    model.parameters.template get<TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)]  = 0.15;
    model.parameters.template get<RelativeTransmissionNoSymptoms<double>>()[mio::AgeGroup(0)]    = 0.5;
    model.parameters.template get<RiskOfInfectionFromSymptomatic<double>>()[mio::AgeGroup(0)]    = 0.0;
    model.parameters.template get<MaxRiskOfInfectionFromSymptomatic<double>>()[mio::AgeGroup(0)] = 0.4;
    model.parameters.template get<RecoveredPerInfectedNoSymptoms<double>>()[mio::AgeGroup(0)]    = 0.2;
    model.parameters.template get<SeverePerInfectedSymptoms<double>>()[mio::AgeGroup(0)]         = 0.1;
    model.parameters.template get<CriticalPerSevere<double>>()[mio::AgeGroup(0)]                 = 0.1;
    model.parameters.template get<DeathsPerCritical<double>>()[mio::AgeGroup(0)]                 = 0.1;

    model.parameters.template get<ReducExposedPartialImmunity<double>>()[mio::AgeGroup(0)]                     = 0.8;
    model.parameters.template get<ReducExposedImprovedImmunity<double>>()[mio::AgeGroup(0)]                    = 0.331;
    model.parameters.template get<ReducInfectedSymptomsPartialImmunity<double>>()[mio::AgeGroup(0)]            = 0.65;
    model.parameters.template get<ReducInfectedSymptomsImprovedImmunity<double>>()[mio::AgeGroup(0)]           = 0.243;
    model.parameters.template get<ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[mio::AgeGroup(0)]  = 0.1;
    model.parameters.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[mio::AgeGroup(0)] = 0.091;
    model.parameters.template get<ReducTimeInfectedMild<double>>()[mio::AgeGroup(0)]                           = 0.9;

    model.parameters.template get<ICUCapacity<double>>()          = 100;
    model.parameters.template get<TestAndTraceCapacity<double>>() = 0.0143;

    const size_t daily_vaccinations = 0;
    const size_t num_days           = 300;
    model.parameters.template get<DailyPartialVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.template get<DailyFullVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.template get<DailyBoosterVaccinations<double>>().resize(mio::SimulationDay(num_days));
    for (size_t i = 0; i < num_days; ++i) {
        auto num_vac = static_cast<double>(i * daily_vaccinations);
        model.parameters.template get<DailyPartialVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(i)}] =
            num_vac;
        model.parameters.template get<DailyFullVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(i)}] =
            num_vac;
        model.parameters.template get<DailyBoosterVaccinations<double>>()[{mio::AgeGroup(0), mio::SimulationDay(i)}] =
            num_vac;
    }

    mio::ContactMatrixGroup& contact_matrix = model.parameters.template get<ContactPatterns<double>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.0));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant(1, 1, 0.7), mio::SimulationTime(30.));
    model.parameters.template get<Seasonality<double>>() = 0.2;
}

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    double t0   = 0;
    double tmax = 2;

    std::vector<double> dts = {2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625};
    const double dt_ref     = 1e-6;

    std::vector<double> base_pop((size_t)InfectionState::Count, 0.0);
    base_pop[(size_t)InfectionState::ExposedNaive]                       = 20;
    base_pop[(size_t)InfectionState::ExposedImprovedImmunity]            = 20;
    base_pop[(size_t)InfectionState::ExposedPartialImmunity]             = 20;
    base_pop[(size_t)InfectionState::InfectedNoSymptomsNaive]            = 30;
    base_pop[(size_t)InfectionState::InfectedNoSymptomsPartialImmunity]  = 30;
    base_pop[(size_t)InfectionState::InfectedNoSymptomsImprovedImmunity] = 30;
    base_pop[(size_t)InfectionState::InfectedSymptomsNaive]              = 40;
    base_pop[(size_t)InfectionState::InfectedSymptomsPartialImmunity]    = 40;
    base_pop[(size_t)InfectionState::InfectedSymptomsImprovedImmunity]   = 40;
    base_pop[(size_t)InfectionState::InfectedSevereNaive]                = 30;
    base_pop[(size_t)InfectionState::InfectedSevereImprovedImmunity]     = 30;
    base_pop[(size_t)InfectionState::InfectedSeverePartialImmunity]      = 30;
    base_pop[(size_t)InfectionState::InfectedCriticalNaive]              = 20;
    base_pop[(size_t)InfectionState::InfectedCriticalPartialImmunity]    = 20;
    base_pop[(size_t)InfectionState::InfectedCriticalImprovedImmunity]   = 20;
    base_pop[(size_t)InfectionState::SusceptibleNaive]                   = 1000;
    base_pop[(size_t)InfectionState::SusceptiblePartialImmunity]         = 1200;
    base_pop[(size_t)InfectionState::SusceptibleImprovedImmunity]        = 1000;
    base_pop[(size_t)InfectionState::TemporaryImmunePartialImmunity]     = 60;
    base_pop[(size_t)InfectionState::TemporaryImmuneImprovedImmunity]    = 70;

    double frac_non_commuter = 0.7;
    double frac_commuter     = 0.3;

    // ========================================================================
    // Reference solution with lagrangian
    // ========================================================================
    std::cout << "Calculating exact Lagrangian Reference (dt = " << dt_ref << ")..." << std::flush;

    StandardModelLagrangianSecirts lag_model(1, 1); // 1 AgeGroup, 1 Commuter Group
    setup_parameters(lag_model);

    for (size_t state = 0; state < (size_t)InfectionState::Count; ++state) {
        lag_model.populations[{mio::AgeGroup(0), CommuterType::NonCommuter, (InfectionState)state}] =
            base_pop[state] * frac_non_commuter;
        lag_model.populations[{mio::AgeGroup(0), CommuterType::CommuterBase, (InfectionState)state}] =
            base_pop[state] * frac_commuter;
    }

    mio::FlowSimulation<double, StandardModelLagrangianSecirts> sim_lag(lag_model, t0, dt_ref);
    auto integrator_ref = std::make_shared<mio::ExplicitStepperWrapper<double, boost::numeric::odeint::runge_kutta4>>();
    sim_lag.set_integrator(integrator_ref);
    sim_lag.advance(tmax);

    Eigen::VectorXd lag_final_full     = sim_lag.get_result().get_last_value();
    Eigen::VectorXd ref_final_commuter = lag_final_full.segment(29, 29);
    std::cout << " Done!\n" << std::endl;

    // ========================================================================
    // 2. Stage-Aligned Setup
    // ========================================================================
    Model<double> sa_model(1);
    setup_parameters(sa_model);

    Eigen::VectorXd initial_totals(29);
    Eigen::VectorXd initial_commuter(29);

    for (size_t state = 0; state < 29; ++state) {
        sa_model.populations[{mio::AgeGroup(0), (InfectionState)state}] = base_pop[state];
        initial_totals[state]                                           = base_pop[state];
        initial_commuter[state]                                         = base_pop[state] * frac_commuter;
    }

    std::vector<std::string> method_names = {"RK1", "RK2", "RK3", "RK4"};

    for (const auto& method_name : method_names) {
        std::cout << "========== Testing Stage-Aligned " << method_name << " against Lagrangian ==========\n";

        double prev_abs_error = -1.0;
        double prev_dt        = -1.0;

        for (auto dt : dts) {
            Eigen::VectorXd current_totals   = initial_totals;
            Eigen::VectorXd current_commuter = initial_commuter;

            for (double t = t0; t < tmax; t += dt) {
                if (t + dt > tmax + 1e-10)
                    break;

                if (method_name == "RK1") {
                    flow_based_mobility_returns_secirts_rk1(current_commuter, current_totals, sa_model, t, dt);
                }
                else if (method_name == "RK2") {
                    flow_based_mobility_returns_secirts_rk2(current_commuter, current_totals, sa_model, t, dt);
                }
                else if (method_name == "RK3") {
                    flow_based_mobility_returns_secirts_rk3(current_commuter, current_totals, sa_model, t, dt);
                }
                else if (method_name == "RK4") {
                    flow_based_mobility_returns_secirts_rk4(current_commuter, current_totals, sa_model, t, dt);
                }
            }

            Eigen::VectorXd error = current_commuter - ref_final_commuter;
            double max_abs_error  = error.cwiseAbs().maxCoeff();

            double max_rel_error = 0.0;
            for (int i = 0; i < 29; ++i) {
                if (std::abs(ref_final_commuter[i]) > 1e-10) {
                    double rel_err = std::abs(error[i]) / std::abs(ref_final_commuter[i]);
                    max_rel_error  = std::max(max_rel_error, rel_err);
                }
            }

            std::string order_str = "   ---";
            if (prev_abs_error > 0.0 && max_abs_error > 0.0 && prev_dt > 0.0) {
                double order = std::log(prev_abs_error / max_abs_error) / std::log(prev_dt / dt);
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(2) << order;
                order_str = oss.str();
            }

            std::cout << "  dt = " << std::scientific << std::setprecision(6) << dt
                      << " | Max abs error: " << max_abs_error << ", Max rel error: " << max_rel_error
                      << " | Order: " << order_str << "\n";

            prev_abs_error = max_abs_error;
            prev_dt        = dt;
        }
        std::cout << std::endl;
    }

    return 0;
}