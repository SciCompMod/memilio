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
#include "lct_secir/initializer_flows.h"

#include "memilio/config.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <iostream>
#include <vector>

namespace params
{
// num_subcompartments is used as a template argument and has to be a constexpr.
constexpr int num_subcompartments = 3;
constexpr size_t num_groups       = 6;

// Parameters
const ScalarType dt                             = 0.01;
const ScalarType seasonality                    = 0.;
const ScalarType RelativeTransmissionNoSymptoms = 1.;
const ScalarType RiskOfInfectionFromSymptomatic = 0.3;
const ScalarType age_group_sizes[]              = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};

const ScalarType TransmissionProbabilityOnContact_age[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

const ScalarType TimeExposed_age[]            = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
const ScalarType TimeInfectedNoSymptoms_age[] = {2.74, 2.74, 2.565, 2.565, 2.565, 2.565};
const ScalarType TimeInfectedSymptoms_age[]   = {7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775};
const ScalarType TimeInfectedSevere_age[]     = {5., 5., 5.925, 7.55, 8.5, 11.};
const ScalarType TimeInfectedCritical_age[]   = {6.95, 6.95, 6.86, 17.36, 17.1, 11.6};

const ScalarType RecoveredPerInfectedNoSymptoms_age[] = {1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8};
const ScalarType SeverePerInfectedSymptoms_age[]      = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType CriticalPerSevere_age[]              = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType DeathsPerCritical_age[]              = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};
} // namespace params

std::vector<ScalarType> get_initial_values(size_t num_subcomp)
{
    const std::vector<ScalarType> init_compartments = {8.28311e+07, 13489.3, 10468,  22297.6,
                                                       1838.56,     571.463, 275292, 0};
    std::vector<ScalarType> initial_value_vector;
    initial_value_vector.push_back(init_compartments[(int)mio::lsecir::InfectionState::Susceptible]);
    // Distribute value equally to the subcompartments.
    for (size_t i = 0; i < num_subcomp; i++) {
        initial_value_vector.push_back(init_compartments[(int)mio::lsecir::InfectionState::Exposed] / num_subcomp);
    }
    for (size_t i = 0; i < num_subcomp; i++) {
        initial_value_vector.push_back(init_compartments[(int)mio::lsecir::InfectionState::InfectedNoSymptoms] /
                                       num_subcomp);
    }
    for (size_t i = 0; i < num_subcomp; i++) {
        initial_value_vector.push_back(init_compartments[(int)mio::lsecir::InfectionState::InfectedSymptoms] /
                                       num_subcomp);
    }
    for (size_t i = 0; i < num_subcomp; i++) {
        initial_value_vector.push_back(init_compartments[(int)mio::lsecir::InfectionState::InfectedSevere] /
                                       num_subcomp);
    }
    for (size_t i = 0; i < num_subcomp; i++) {
        initial_value_vector.push_back(init_compartments[(int)mio::lsecir::InfectionState::InfectedCritical] /
                                       num_subcomp);
    }
    initial_value_vector.push_back(init_compartments[(int)mio::lsecir::InfectionState::Recovered]);
    initial_value_vector.push_back(init_compartments[(int)mio::lsecir::InfectionState::Dead]);

    return initial_value_vector;
}

/** 
* @brief Perform a fictive simulation with realistic parameters and contacts, such that the reproduction number 
*   is approximately 1 at the beginning and rising or dropping at simulation time 2.
*   
*   This scenario should enable a comparison of the qualitative behavior of different LCT models.
*   
* @param[in] R0 Define R0 from simulation time 2 on. Please use a number > 0.
* @param[in] tmax End time of the simulation.
* @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
* @returns Any io errors that happen during saving the results.
*/
mio::IOResult<void> simulate_lct_model(ScalarType R0, ScalarType tmax, bool save_subcompartments,
                                       std::string save_dir = "")
{
    using namespace params;
    std::cout << "Simulation with LCT model and " << num_subcompartments << " subcompartments." << std::endl;

    // Initialize model.
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
    ScalarType total_population = 0;
    for (size_t group = 0; group < num_groups; group++) {
        total_population += age_group_sizes[group];
    }
    ScalarType TimeExposed                      = 0;
    ScalarType TimeInfectedNoSymptoms           = 0;
    ScalarType TimeInfectedSymptoms             = 0;
    ScalarType TimeInfectedSevere               = 0;
    ScalarType TimeInfectedCritical             = 0;
    ScalarType TransmissionProbabilityOnContact = 0;
    ScalarType RecoveredPerInfectedNoSymptoms   = 0;
    ScalarType SeverePerInfectedSymptoms        = 0;
    ScalarType CriticalPerSevere                = 0;
    ScalarType DeathsPerCritical                = 0;
    for (size_t group = 0; group < num_groups; group++) {
        TimeExposed += age_group_sizes[group] * TimeExposed_age[group] / total_population;
        TimeInfectedNoSymptoms += age_group_sizes[group] * TimeInfectedNoSymptoms_age[group] / total_population;
        TimeInfectedSymptoms += age_group_sizes[group] * TimeInfectedSymptoms_age[group] / total_population;
        TimeInfectedSevere += age_group_sizes[group] * TimeInfectedSevere_age[group] / total_population;
        TimeInfectedCritical += age_group_sizes[group] * TimeInfectedCritical_age[group] / total_population;
        TransmissionProbabilityOnContact +=
            age_group_sizes[group] * TransmissionProbabilityOnContact_age[group] / total_population;
        RecoveredPerInfectedNoSymptoms +=
            age_group_sizes[group] * RecoveredPerInfectedNoSymptoms_age[group] / total_population;
        SeverePerInfectedSymptoms += age_group_sizes[group] * SeverePerInfectedSymptoms_age[group] / total_population;
        CriticalPerSevere += age_group_sizes[group] * CriticalPerSevere_age[group] / total_population;
        DeathsPerCritical += age_group_sizes[group] * DeathsPerCritical_age[group] / total_population;
    }

    model.parameters.get<mio::lsecir::TimeExposed>()[0]                      = TimeExposed;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0]           = TimeInfectedNoSymptoms;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]             = TimeInfectedSymptoms;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]               = TimeInfectedSevere;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]             = TimeInfectedCritical;
    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = TransmissionProbabilityOnContact;

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = RelativeTransmissionNoSymptoms;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = RiskOfInfectionFromSymptomatic;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = RecoveredPerInfectedNoSymptoms;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = SeverePerInfectedSymptoms;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = CriticalPerSevere;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = DeathsPerCritical;

    ScalarType contacts_R1 =
        1. / (TransmissionProbabilityOnContact *
              (TimeInfectedNoSymptoms * RelativeTransmissionNoSymptoms +
               (1 - RecoveredPerInfectedNoSymptoms) * TimeInfectedSymptoms * RiskOfInfectionFromSymptomatic));

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    if (R0 <= 1.) {
        // Perform simulation with dropping R0.
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, contacts_R1));
        contact_matrix[0].add_damping(0., mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(R0, mio::SimulationTime(2.));
    }
    else {
        // Perform simulation with rising R0.
        contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, R0 * contacts_R1));
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(-1.));
        contact_matrix[0].add_damping(1 - 1. / R0, mio::SimulationTime(1.9));
        contact_matrix[0].add_damping(0., mio::SimulationTime(2.));
    }

    model.parameters.get<mio::lsecir::ContactPatterns>() = mio::UncertainContactMatrix<ScalarType>(contact_matrix);
    model.parameters.get<mio::lsecir::Seasonality>()     = seasonality;

    auto initial_values = get_initial_values(num_subcompartments);
    for (size_t i = 0; i < model.populations.get_num_compartments(); i++) {
        model.populations[i] = initial_values[i];
    }

    // Perform simulation.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max so that we have a fixed time step and can compare to the result with one group.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
    std::cout << result[0] << std::endl;
    // Calculate result.
    mio::TimeSeries<ScalarType> populations = model.calculate_compartments(result);

    if (!save_dir.empty()) {

        std::string R0string = std::to_string(R0);
        std::string filename = save_dir + "fictional_lct_" + R0string.substr(0, R0string.find(".") + 2) + "_" +
                               std::to_string(num_subcompartments);
        if (save_subcompartments) {
            filename                               = filename + "_subcompartments.h5";
            auto result_interpolated               = mio::interpolate_simulation_result(result);
            mio::IOResult<void> save_result_status = mio::save_result({result_interpolated}, {0}, 1, filename);
        }
        else {
            filename                               = filename + ".h5";
            mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename);
        }
    }
    // std::cout << "Final size: " << std::fixed << std::setprecision(6)
    //           << total_population - populations.get_last_value()[0] << std::endl;
    // std::cout << std::endl;
    // order: number of subcompartments: 1,3,10,50
    //const ScalarType erg[] = {66187880.970838, 66177693.857084, 66173548.328663, 66172040.662903}; //R=2
    //const ScalarType erg[] = {81489438.000375, 81487771.513275, 81487273.311926, 81487137.403808}; //R=4
    // const ScalarType erg[] = {83151138.102138, 83151130.435465, 83151128.866535, 83151128.512536}; //R=10
    // std::cout << "Absolute deviation: " << std::endl;
    // for (int i = 0; i < 4; i++) {
    //     std::cout << "i= " << i << ": " << (erg[i] - erg[0]) << std::endl;
    // }
    // std::cout << "Relative deviation: " << std::endl;
    // for (int i = 0; i < 4; i++) {
    //     std::cout << "i= " << i << ": " << (erg[i] - erg[0]) / erg[0] << std::endl;
    // }
    return mio::success();
}

int main(int argc, char** argv)
{
    std::string save_dir       = "../../data/simulation_lct_noage/riseR0short/";
    ScalarType R0              = 2.;
    bool save_subcompartments  = false;
    ScalarType simulation_days = 12;
    if (argc > 2) {
        R0              = std::stod(argv[1]);
        simulation_days = std::stod(argv[2]);
    }
    if (argc > 3) {
        save_dir = argv[3];
    }
    auto result = simulate_lct_model(R0, simulation_days, save_subcompartments, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
