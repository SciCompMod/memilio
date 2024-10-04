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
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/io/mobility_io.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <iostream>
#include <vector>
namespace params
{
// Necessary because num_subcompartments is used as a template argument and has to be a constexpr.
constexpr int num_subcompartments = 10;
constexpr size_t num_groups       = 6;

// Parameters
const ScalarType dt                             = 0.01;
const ScalarType seasonality                    = 0.;
const ScalarType RelativeTransmissionNoSymptoms = 1.;
const ScalarType RiskOfInfectionFromSymptomatic = 0.3;
const ScalarType age_group_sizes[]              = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
const ScalarType total_population               = 83155031.0;

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

mio::TimeSeries<ScalarType> add_age_groups(mio::TimeSeries<ScalarType> ageres)
{
    using namespace params;
    size_t infstatecount = (size_t)mio::lsecir::InfectionState::Count;
    mio::TimeSeries<ScalarType> noage(infstatecount);

    for (Eigen::Index timepoint = 0; timepoint < ageres.get_num_time_points(); ++timepoint) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(infstatecount);
        for (size_t infstate = 0; infstate < infstatecount; infstate++) {
            for (size_t group = 0; group < num_groups; group++) {
                result[infstate] += ageres.get_value(timepoint)[group * infstatecount + infstate];
            }
        }
        noage.add_time_point(ageres.get_time(timepoint), result);
    }

    return noage;
}

/** @brief Returns contact matrix in relation to defined R0.
* Contacts are defined such that R0 equals 1 at the beginning of the simulation and jumps to R0 in 
* the time interval [1.9,2.0].
*/
mio::IOResult<mio::UncertainContactMatrix<ScalarType>> get_contact_matrix(bool age = true)
{
    using namespace params;
    std::string contact_data_dir          = "../../data/contacts/";
    const std::string contact_locations[] = {"home", "school_pf_eig", "work", "other"};
    const size_t num_locations            = 4;
    size_t mat_size                       = num_groups;
    if (!age) {
        mat_size = 1;
    }
    auto contact_matrices = mio::ContactMatrixGroup(num_locations, mat_size);
    // Load and set baseline contacts for each contact location.
    for (size_t location = 0; location < num_locations; location++) {
        BOOST_OUTCOME_TRY(auto&& baseline, mio::read_mobility_plain(contact_data_dir + "baseline_" +
                                                                    contact_locations[location] + ".txt"));
        if (!age) {
            ScalarType base = 0.;
            for (size_t i = 0; i < num_groups; i++) {
                for (size_t j = 0; j < num_groups; j++) {
                    // Calculate a weighted average according to the age group sizes of the total contacts.
                    base += age_group_sizes[i] / total_population * baseline(i, j);
                }
            }
            contact_matrices[location].get_baseline() = Eigen::MatrixXd::Constant(1, 1, base);
            contact_matrices[location].get_minimum()  = Eigen::MatrixXd::Zero(1, 1);
        }
        else {
            contact_matrices[location].get_baseline() = baseline;
            contact_matrices[location].get_minimum()  = Eigen::MatrixXd::Zero(num_groups, num_groups);
        }
    }
    std::cout << "Contact Pattern used:" << std::endl;
    std::cout << mio::UncertainContactMatrix<ScalarType>(contact_matrices).get_cont_freq_mat().get_matrix_at(0)
              << std::endl;
    std::cout << std::endl;
    return mio::success(mio::UncertainContactMatrix<ScalarType>(contact_matrices));
}

std::vector<std::vector<ScalarType>> get_initialization(bool age = false, size_t agegroup_init = 0)
{
    using namespace params;
    ScalarType num_exposed = 100.;
    std::vector<std::vector<ScalarType>> init;
    if (age) {
        for (size_t group = 0; group < num_groups; group++) {
            std::vector<ScalarType> init_vector_group({total_population, 0., 0., 0., 0., 0., 0., 0.});
            if (group == agegroup_init) {
                init_vector_group[0] = total_population - num_exposed;
                init_vector_group[1] = num_exposed;
            }
            init.push_back(init_vector_group);
        }
    }
    else {
        init.push_back({total_population - num_exposed, num_exposed, 0., 0., 0., 0., 0., 0.});
    }
    return init;
}

/** 
* @brief Perform a fictive simulation with realistic parameters and contacts, such that the reproduction number 
*   is approximately 1 at the beginning and rising or dropping at simulation time 2.
*   
*   This scenario should enable a comparison of the qualitative behavior of different LCT models.
*   
* @param[in] tmax End time of the simulation.
* @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
* @returns Any io errors that happen during saving the results.
*/
mio::IOResult<void> simulate_ageres_model(size_t agegroup_init, ScalarType tmax, std::string save_dir = "")
{
    using namespace params;
    std::cout << "Simulation with " << num_subcompartments
              << " subcompartments for the initial exposed population in age group " << agegroup_init << "."
              << std::endl;
    // ----- Initialize age resolved model. -----
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState, LctState, LctState, LctState, LctState, LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
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
    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix());
    model.parameters.get<mio::lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.get<mio::lsecir::Seasonality>()     = seasonality;

    auto init_vector = get_initialization(true, agegroup_init);
    for (size_t group = 0; group < num_groups; group++) {
        model.populations[group * LctState::Count + 0]                   = init_vector[group][0]; //S
        model.populations[group * LctState::Count + LctState::Count - 2] = init_vector[group][6]; //R
        model.populations[group * LctState::Count + LctState::Count - 1] = init_vector[group][7]; //D
        for (size_t i = 1; i < (size_t)InfState::Count - 2; i++) {
            for (size_t subcomp = 0; subcomp < num_subcompartments; subcomp++) {
                model.populations[group * LctState::Count + (i - 1) * num_subcompartments + 1 + subcomp] =
                    init_vector[group][i] / num_subcompartments;
            }
        }
    }

    // Perform simulation.
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(
        0, tmax, dt, model,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            1e-10, 1e-5, 0, dt));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations = add_age_groups(model.calculate_compartments(result));

    if (!save_dir.empty()) {
        std::string filename = save_dir + "fictional_lct_ageres_" + std::to_string(num_subcompartments) +
                               "_agegroupinit_" + std::to_string(agegroup_init) + ".h5";
        mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename);
    }

    return mio::success();
}

mio::IOResult<void> simulate_notageres_model(ScalarType tmax, std::string save_dir = "")
{
    using namespace params;
    std::cout << "Simulation with " << num_subcompartments << " subcompartments without age resolution." << std::endl;
    // ----- Initialize age resolved model. -----
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
    ScalarType TimeExposed_noage                      = 0;
    ScalarType TimeInfectedNoSymptoms_noage           = 0;
    ScalarType TimeInfectedSymptoms_noage             = 0;
    ScalarType TimeInfectedSevere_noage               = 0;
    ScalarType TimeInfectedCritical_noage             = 0;
    ScalarType TransmissionProbabilityOnContact_noage = 0;
    ScalarType RecoveredPerInfectedNoSymptoms_noage   = 0;
    ScalarType SeverePerInfectedSymptoms_noage        = 0;
    ScalarType CriticalPerSevere_noage                = 0;
    ScalarType DeathsPerCritical_noage                = 0;
    for (size_t group = 0; group < num_groups; group++) {
        TimeExposed_noage += age_group_sizes[group] * TimeExposed[group] / total_population;
        TimeInfectedNoSymptoms_noage += age_group_sizes[group] * TimeInfectedNoSymptoms[group] / total_population;
        TimeInfectedSymptoms_noage += age_group_sizes[group] * TimeInfectedSymptoms[group] / total_population;
        TimeInfectedSevere_noage += age_group_sizes[group] * TimeInfectedSevere[group] / total_population;
        TimeInfectedCritical_noage += age_group_sizes[group] * TimeInfectedCritical[group] / total_population;
        TransmissionProbabilityOnContact_noage +=
            age_group_sizes[group] * TransmissionProbabilityOnContact[group] / total_population;
        RecoveredPerInfectedNoSymptoms_noage +=
            age_group_sizes[group] * RecoveredPerInfectedNoSymptoms[group] / total_population;
        SeverePerInfectedSymptoms_noage += age_group_sizes[group] * SeverePerInfectedSymptoms[group] / total_population;
        CriticalPerSevere_noage += age_group_sizes[group] * CriticalPerSevere[group] / total_population;
        DeathsPerCritical_noage += age_group_sizes[group] * DeathsPerCritical[group] / total_population;
    }

    model.parameters.get<mio::lsecir::TimeExposed>()[0]                      = TimeExposed_noage;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0]           = TimeInfectedNoSymptoms_noage;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]             = TimeInfectedSymptoms_noage;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]               = TimeInfectedSevere_noage;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]             = TimeInfectedCritical_noage;
    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = TransmissionProbabilityOnContact_noage;

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = RelativeTransmissionNoSymptoms;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = RiskOfInfectionFromSymptomatic;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = RecoveredPerInfectedNoSymptoms_noage;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = SeverePerInfectedSymptoms_noage;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = CriticalPerSevere_noage;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = DeathsPerCritical_noage;

    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(false));
    model.parameters.get<mio::lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.get<mio::lsecir::Seasonality>()     = seasonality;

    auto init_vector                       = get_initialization(false);
    model.populations[0]                   = init_vector[0][0]; //S
    model.populations[LctState::Count - 2] = init_vector[0][6]; //R
    model.populations[LctState::Count - 1] = init_vector[0][7]; //D
    for (size_t i = 1; i < (size_t)InfState::Count - 2; i++) {
        for (size_t subcomp = 0; subcomp < num_subcompartments; subcomp++) {
            model.populations[(i - 1) * num_subcompartments + 1 + subcomp] = init_vector[0][i] / num_subcompartments;
        }
    }

    // Perform simulation.
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(
        0, tmax, dt, model,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>(
            1e-10, 1e-5, 0, dt));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations = model.calculate_compartments(result);

    if (!save_dir.empty()) {
        std::string filename = save_dir + "fictional_lct_notageres_" + std::to_string(num_subcompartments) + ".h5";
        mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename);
    }

    return mio::success();
}

int main()
{
    std::string save_dir = "../../data/simulation_lct_agevsnoage/";
    ScalarType tmax      = 40;
    // Simulation with initial exposed population in age group 1.
    auto result = simulate_ageres_model(1, tmax, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    // Simulation with initial exposed population in age group 5.
    result = simulate_ageres_model(5, tmax, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    // Simulation without age resolution.
    result = simulate_notageres_model(tmax, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
