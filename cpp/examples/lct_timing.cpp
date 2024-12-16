/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Lena Ploetzke
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
#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"

#include "memilio/config.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <iostream>
#include <vector>
#include <omp.h>

namespace params
{
// Define epidemiological parameters and parameters needed for the simulation.
constexpr size_t num_subcompartments = NUM_SUBCOMPARTMENTS;
constexpr size_t num_groups          = 1;

const ScalarType dt                = 0.01;
const ScalarType age_group_sizes[] = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
const ScalarType total_population  = 83155031.0;

const ScalarType seasonality                        = 0.;
const ScalarType RelativeTransmissionNoSymptoms     = 1.;
const ScalarType RiskOfInfectionFromSymptomatic     = 0.3;
const ScalarType TransmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

const ScalarType TimeExposed[]            = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
const ScalarType TimeInfectedNoSymptoms[] = {2.74, 2.74, 2.565, 2.565, 2.565, 2.565};
const ScalarType TimeInfectedSymptoms[]   = {7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775};
const ScalarType TimeInfectedSevere[]     = {5., 5., 5.925, 7.55, 8.5, 11.};
const ScalarType TimeInfectedCritical[]   = {6.95, 6.95, 6.86, 17.36, 17.1, 11.6};

const ScalarType RecoveredPerInfectedNoSymptoms[] = {1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8};
const ScalarType SeverePerInfectedSymptoms[]      = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType CriticalPerSevere[]              = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType DeathsPerCritical[]              = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};

} // namespace params

/** @brief Initial value vector for the simulation.
*   It is assumed that all age groups use equal LctStates.
* @tparam LctStates LctState of all the age groups.
*/
std::vector<ScalarType> get_initial_values()
{
    using namespace params;
    // Vector is a "random vector" taken from another example. Just need some realistic values.
    const std::vector<std::vector<ScalarType>> init_compartments = {
        {3966564.2110, 664.2367, 545.5523, 1050.3946, 5.6045, 0.5844, 307.4165, 0.},
        {7500988.3044, 2108.8502, 1732.0453, 3334.8427, 17.7934, 1.8555, 478.3085, 0.},
        {18874457.8051, 12584.3371, 9674.4579, 21348.6877, 340.0557, 29.5323, 2857.1243, 0.},
        {28612752.3265, 13788.4953, 10600.1783, 22967.4792, 1537.0744, 530.3313, 3990.1151, 0.},
        {8134534.1310, 4612.1712, 3545.6978, 7567.8082, 1553.0877, 937.6034, 588.5007, 0.},
        {928318.1474, 1595.6020, 1226.6506, 2595.1188, 948.2972, 400.0181, 1350.1658, 0.}};
    std::vector<ScalarType> initial_value_vector;
    for (size_t age = 0; age < num_groups; age++) {
        std::vector<ScalarType> init_age;
        init_age.push_back(init_compartments[age][(int)mio::lsecir::InfectionState::Susceptible]);
        // Distribute value equally to the subcompartments.
        for (size_t i = 0; i < num_subcompartments; i++) {
            init_age.push_back(init_compartments[age][(int)mio::lsecir::InfectionState::Exposed] / num_subcompartments);
        }
        for (size_t i = 0; i < num_subcompartments; i++) {
            init_age.push_back(init_compartments[age][(int)mio::lsecir::InfectionState::InfectedNoSymptoms] /
                               num_subcompartments);
        }
        for (size_t i = 0; i < num_subcompartments; i++) {
            init_age.push_back(init_compartments[age][(int)mio::lsecir::InfectionState::InfectedSymptoms] /
                               num_subcompartments);
        }
        for (size_t i = 0; i < num_subcompartments; i++) {
            init_age.push_back(init_compartments[age][(int)mio::lsecir::InfectionState::InfectedSevere] /
                               num_subcompartments);
        }
        for (size_t i = 0; i < num_subcompartments; i++) {
            init_age.push_back(init_compartments[age][(int)mio::lsecir::InfectionState::InfectedCritical] /
                               num_subcompartments);
        }
        init_age.push_back(init_compartments[age][(int)mio::lsecir::InfectionState::Recovered]);
        init_age.push_back(init_compartments[age][(int)mio::lsecir::InfectionState::Dead]);

        initial_value_vector.insert(initial_value_vector.end(), init_age.begin(), init_age.end());
    }
    return initial_value_vector;
}

/**
 * @brief Performs multiple simulations with one model to get an average run time.
 * @tparam num_subcompartments number of subcompartments used for all compartments and all age groups.
 */
void simulate(size_t num_warm_up_runs, size_t num_runs, ScalarType tmax)
{
    using namespace params;
    std::cout << "{ \"Agegroups\": " << num_groups << ",\n\"Subcompartments\": " << num_subcompartments << ", "
              << std::endl;
    // ----- Initialize age resolved model. -----
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Define epidemiological parameters.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[group]            = TimeExposed[group];
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[group] = TimeInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[group]   = TimeInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[group]     = TimeInfectedSevere[group];
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[group]   = TimeInfectedCritical[group];
        model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[group] =
            TransmissionProbabilityOnContact[group];

        model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[group] = RelativeTransmissionNoSymptoms;
        model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[group] = RiskOfInfectionFromSymptomatic;

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[group] =
            RecoveredPerInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[group] = SeverePerInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[group]         = CriticalPerSevere[group];
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[group]         = DeathsPerCritical[group];
    }
    // Realistic contacts.
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::lsecir::ContactPatterns>();
    Eigen::MatrixXd contact_matrix_eigen(6, 6);
    contact_matrix_eigen << 3.9547, 1.1002, 2.9472, 2.05, 0.3733, 0.0445, 0.3327, 3.5892, 1.236, 1.9208, 0.2681, 0.0161,
        0.246, 0.7124, 5.6518, 3.2939, 0.2043, 0.0109, 0.1742, 0.8897, 3.3124, 4.5406, 0.4262, 0.0214, 0.0458, 0.1939,
        0.5782, 1.3825, 1.473, 0.0704, 0.1083, 0.1448, 0.4728, 0.9767, 0.6266, 0.1724;
    contact_matrix[0] = mio::ContactMatrix(contact_matrix_eigen.block(0, 0, (size_t)num_groups, (size_t)num_groups));

    model.parameters.get<mio::lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.get<mio::lsecir::Seasonality>()     = seasonality;

    // Set initial values;
    auto initial_values = get_initial_values();
    for (size_t i = 0; i < model.populations.get_num_compartments(); i++) {
        model.populations[i] = initial_values[i];
    }
    // Integrator.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max so that we have a fixed time step.
    // integrator->set_dt_min(dt);
    // integrator->set_dt_max(dt);

    // Warm up runs.
    mio::set_log_level(mio::LogLevel::off);
    for (size_t i = 0; i < num_warm_up_runs; i++) {
        mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
    }
    auto result = mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
    std::cout << "\"Steps\": " << result.get_num_time_points() << "," << std::endl;

    // Runs with timing.
    ScalarType total = 0;
    for (size_t i = 0; i < num_runs; i++) {
        total -= omp_get_wtime();
        mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
        total += omp_get_wtime();
    }
    std::cout << "\"Time\": " << total / num_runs << "\n}," << std::endl;
    mio::set_log_level(mio::LogLevel::warn);
}

int main(int argc, char** argv)
{
    const ScalarType tmax = 20;
    size_t warm_up        = 10;
    size_t num_runs       = 100;
    if (argc > 2) {
        warm_up  = std::stod(argv[1]);
        num_runs = std::stod(argv[2]);
    }
    simulate(warm_up, num_runs, tmax);
    return 0;
}
