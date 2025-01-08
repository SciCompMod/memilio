/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/math/stepper_wrapper.h"
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
// num_subcompartments is used as a template argument and has to be a constexpr.
constexpr int num_subcompartments = 10;
constexpr size_t num_groups       = 6;

// Define (age-resolved) parameters.
const ScalarType dt                             = 0.01;
const ScalarType seasonality                    = 0.;
const ScalarType relativeTransmissionNoSymptoms = 1.;
const ScalarType riskOfInfectionFromSymptomatic = 0.3;
const ScalarType age_group_sizes[]              = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
const ScalarType total_population               = 83155031.0;

const ScalarType transmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

const ScalarType timeExposed[]            = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
const ScalarType timeInfectedNoSymptoms[] = {2.74, 2.74, 2.565, 2.565, 2.565, 2.565};
const ScalarType timeInfectedSymptoms[]   = {7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775};
const ScalarType timeInfectedSevere[]     = {5., 5., 5.925, 7.55, 8.5, 11.};
const ScalarType timeInfectedCritical[]   = {6.95, 6.95, 6.86, 17.36, 17.1, 11.6};

const ScalarType recoveredPerInfectedNoSymptoms[] = {1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8};
const ScalarType severePerInfectedSymptoms[]      = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType criticalPerSevere[]              = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType deathsPerCritical[]              = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};
} // namespace params

/** 
* @brief Function to transform an age-resolved simulation result into a result without age resolution.
*
* Sums up the values in the age groups to transform the simulation result into a result without age resolution. 
* For the sake of comparability, we use non-age-resolved results for visualizations.
* This implementation is only valid if the simulation is run with equal LctStates for all groups or if the result 
* does not contain any subcompartments.
*   
* @param[in] ageresolved_result TimeSeries with an age-resolved simulation result.
* @returns TimeSeries with the result where the values of the age groups are summed up.
*/
mio::TimeSeries<ScalarType> sum_age_groups(const mio::TimeSeries<ScalarType> ageresolved_result)
{
    using namespace params;
    size_t infstatecount = size_t((ScalarType)ageresolved_result.get_num_elements() / (ScalarType)num_groups);
    mio::TimeSeries<ScalarType> nonageresolved_result(infstatecount);

    // For each time point, calculate the result without age resolution and add the time point
    // to the non-age-resolved result.
    for (Eigen::Index timepoint = 0; timepoint < ageresolved_result.get_num_time_points(); ++timepoint) {
        Eigen::VectorX<ScalarType> result = Eigen::VectorX<ScalarType>::Zero(infstatecount);
        for (size_t infstate = 0; infstate < infstatecount; infstate++) {
            for (size_t group = 0; group < num_groups; group++) {
                result[infstate] += ageresolved_result.get_value(timepoint)[group * infstatecount + infstate];
            }
        }
        nonageresolved_result.add_time_point(ageresolved_result.get_time(timepoint), result);
    }

    return nonageresolved_result;
}

/** 
* @brief Gets a contact matrix from data files and computes a weighted average for non-age-resolved simulations.
* @param[in] contact_data_dir Directory to the contact data.
* @param[in] resolve_by_age If true, the function gives an age-resolved contact matrix. If false, a weighted average
*    is calculated from the age-resolved data. Default is false.
* @returns The contact matrix or any IO errors that occur during reading the contact data files.
*/
mio::IOResult<mio::UncertainContactMatrix<ScalarType>> get_contact_matrix(std::string contact_data_dir,
                                                                          bool resolve_by_age = false)
{
    using namespace params;
    const std::string contact_locations[] = {"home", "school_pf_eig", "work", "other"};
    const size_t num_locations            = 4;
    size_t matrix_size                    = num_groups;
    if (!resolve_by_age) {
        matrix_size = 1;
    }
    auto contact_matrices = mio::ContactMatrixGroup(num_locations, matrix_size);
    // Load and set baseline contacts for each contact location.
    for (size_t location = 0; location < num_locations; location++) {
        BOOST_OUTCOME_TRY(auto&& baseline, mio::read_mobility_plain(contact_data_dir + "baseline_" +
                                                                    contact_locations[location] + ".txt"));
        if (!resolve_by_age) {
            ScalarType average = 0.;
            for (size_t i = 0; i < num_groups; i++) {
                for (size_t j = 0; j < num_groups; j++) {
                    // Calculate a weighted average according to the age group sizes.
                    average += age_group_sizes[i] / total_population * baseline(i, j);
                }
            }
            contact_matrices[location].get_baseline() = Eigen::MatrixXd::Constant(matrix_size, matrix_size, average);
            contact_matrices[location].get_minimum()  = Eigen::MatrixXd::Zero(matrix_size, matrix_size);
        }
        else {
            contact_matrices[location].get_baseline() = baseline;
            contact_matrices[location].get_minimum()  = Eigen::MatrixXd::Zero(matrix_size, matrix_size);
        }
    }
    return mio::success(mio::UncertainContactMatrix<ScalarType>(contact_matrices));
}

/** 
* @brief Constructs an initial value vector with 100 Exposed individuals and remaining population in Susceptible.
*   
* The function constructs a vector of vectors with one vector for each age group. For each age group, 
* the vector contains initial values for each compartment (without resolution in subcompartments).
* If resolve_by_age is true, the initial number of Exposed individuals in age group agegroup_exposed is set to 100. 
* If resolve_by_age is false, an initial value vector without age resolution is created with 100 Exposed individuals.
* The idea is that all settings for agegroup_exposed are translated into the same setting if we use a model 
*   without age resolution. 
* The vector of vectors should be used as a basis to construct a valid initial value vector with subcompartments. 
*
* @param[in] resolve_by_age If true, the function provides age-resolved data that can be used for initialization 
*   purpose. If false, only one group is used. Default is false.
* @param[in] agegroup_exposed The age group with the 100 initially exposed individuals. 
*    This only makes sense if resolved_by_age is true.
* @returns The initial value vector.
*/
std::vector<std::vector<ScalarType>> get_initialization(bool resolve_by_age = false, size_t agegroup_exposed = 0)
{
    using namespace params;
    ScalarType num_exposed = 100.;
    std::vector<std::vector<ScalarType>> init;
    if (resolve_by_age) {
        for (size_t group = 0; group < num_groups; group++) {
            std::vector<ScalarType> init_vector_group({age_group_sizes[group], 0., 0., 0., 0., 0., 0., 0.});
            if (group == agegroup_exposed) {
                init_vector_group[0] = age_group_sizes[group] - num_exposed;
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
* @brief Perform simulation to examine the impact of the age resolution. 
*
*   The simulation uses LCT models with Covid-19 inspired parameters and a contact rate for Germany. 
*   The simulation is performed with 100 initially Exposed individuals in age group agegroup_exposed. 
*   The remaining compartment sizes are set to zero and the number of Susceptibles is set accordingly.
*   The idea is that all settings for agegroup_exposed are translated into the same setting if we use a model 
*   without age resolution. 
*   Therefore, by comparing different settings for agegroup_exposed, we can assess the impact of age resolution.
*  
* @param[in] agegroup_exposed The age group with the 100 initially exposed individuals. 
* @param[in] tmax End time of the simulation.
* @param[in] contact_data_dir Directory to the contact data.
* @param[in] save_dir Specifies the directory where the results should be stored. 
*   Provide an empty string if the results should not be saved.
* @returns Any io errors that happen during saving the results.
*/
mio::IOResult<void> simulation_with_ageresolution(size_t agegroup_exposed, ScalarType tmax,
                                                  std::string contact_data_dir, std::string save_dir = "")
{
    using namespace params;
    std::cout << "Simulation with " << num_subcompartments
              << " subcompartments and 100 Exposed individuals in age group " << agegroup_exposed << "." << std::endl;
    // Initialize age-resolved model.
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState, LctState, LctState, LctState, LctState, LctState>;
    Model model;

    // Define parameters.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[group]            = timeExposed[group];
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[group] = timeInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[group]   = timeInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[group]     = timeInfectedSevere[group];
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[group]   = timeInfectedCritical[group];
        model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[group] =
            transmissionProbabilityOnContact[group];

        model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[group] = relativeTransmissionNoSymptoms;
        model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[group] = riskOfInfectionFromSymptomatic;

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[group] =
            recoveredPerInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[group] = severePerInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[group]         = criticalPerSevere[group];
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[group]         = deathsPerCritical[group];
    }
    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(contact_data_dir, true));
    model.parameters.get<mio::lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.get<mio::lsecir::Seasonality>()     = seasonality;

    auto init = get_initialization(true, agegroup_exposed);
    // Use init as a basis to define appropriate initial values.
    // Compartment values are distributed equally to subcompartments.
    for (size_t group = 0; group < num_groups; group++) {
        model.populations[group * LctState::Count + 0]                   = init[group][0]; // Susceptible
        model.populations[group * LctState::Count + LctState::Count - 2] = init[group][6]; // Recovered
        model.populations[group * LctState::Count + LctState::Count - 1] = init[group][7]; // Dead
        for (size_t i = 1; i < (size_t)InfState::Count - 2; i++) {
            for (size_t subcomp = 0; subcomp < num_subcompartments; subcomp++) {
                model.populations[group * LctState::Count + (i - 1) * num_subcompartments + 1 + subcomp] =
                    init[group][i] / (ScalarType)num_subcompartments;
            }
        }
    }

    // Set integrator of fifth order with fixed step size and perform simulation.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max to get a fixed step size.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
    // Calculate result without division in subcompartments and without division in age groups.
    mio::TimeSeries<ScalarType> populations = sum_age_groups(model.calculate_compartments(result));

    if (!save_dir.empty()) {
        std::string filename = save_dir + "lct_ageresolved_subcomp" + std::to_string(num_subcompartments) +
                               "_agegroupinit" + std::to_string(agegroup_exposed) + ".h5";
        mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename);
    }

    return mio::success();
}

/** 
* @brief Perform simulation to examine the impact of age resolution without an age resolution for comparison. 
*
*   The simulation uses LCT models with Covid-19 inspired parameters and a contact rate for Germany. 
*   The simulation is performed with 100 initially Exposed individuals. 
*   The remaining compartment sizes are set to zero and the number of Susceptibles is set accordingly.
*   The idea is that all settings of the function simulation_with_ageresolution() are translated into this setting
*   if only a model without age resolution is provided.
*   Therefore, by comparing this result without age resolution with different settings for the function 
*   simulation_with_ageresolution(), we can assess the impact of age resolution.
*  
* @param[in] tmax End time of the simulation.
* @param[in] contact_data_dir Directory to the contact data.
* @param[in] save_dir Specifies the directory where the results should be stored. 
*   Provide an empty string if the results should not be saved.
* @returns Any io errors that happen during saving the results.
*/
mio::IOResult<void> simulation_without_ageresolution(ScalarType tmax, std::string contact_data_dir,
                                                     std::string save_dir = "")
{
    using namespace params;
    std::cout << "Simulation with " << num_subcompartments << " subcompartments without age resolution." << std::endl;
    // Initialize non-age-resolved model (i.e., with one single group).
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Define non-age-resolved parameters by computing weighted averages of the age-resolved parameters.
    ScalarType timeExposed_noage = 0, timeInfectedNoSymptoms_noage = 0, timeInfectedSymptoms_noage = 0,
               timeInfectedSevere_noage = 0, timeInfectedCritical_noage = 0, transmissionProbabilityOnContact_noage = 0,
               recoveredPerInfectedNoSymptoms_noage = 0, severePerInfectedSymptoms_noage = 0,
               criticalPerSevere_noage = 0, deathsPerCritical_noage = 0;
    for (size_t group = 0; group < num_groups; group++) {
        timeExposed_noage += age_group_sizes[group] * timeExposed[group] / total_population;
        timeInfectedNoSymptoms_noage += age_group_sizes[group] * timeInfectedNoSymptoms[group] / total_population;
        timeInfectedSymptoms_noage += age_group_sizes[group] * timeInfectedSymptoms[group] / total_population;
        timeInfectedSevere_noage += age_group_sizes[group] * timeInfectedSevere[group] / total_population;
        timeInfectedCritical_noage += age_group_sizes[group] * timeInfectedCritical[group] / total_population;
        transmissionProbabilityOnContact_noage +=
            age_group_sizes[group] * transmissionProbabilityOnContact[group] / total_population;
        recoveredPerInfectedNoSymptoms_noage +=
            age_group_sizes[group] * recoveredPerInfectedNoSymptoms[group] / total_population;
        severePerInfectedSymptoms_noage += age_group_sizes[group] * severePerInfectedSymptoms[group] / total_population;
        criticalPerSevere_noage += age_group_sizes[group] * criticalPerSevere[group] / total_population;
        deathsPerCritical_noage += age_group_sizes[group] * deathsPerCritical[group] / total_population;
    }

    model.parameters.get<mio::lsecir::TimeExposed>()[0]                      = timeExposed_noage;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0]           = timeInfectedNoSymptoms_noage;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]             = timeInfectedSymptoms_noage;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]               = timeInfectedSevere_noage;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]             = timeInfectedCritical_noage;
    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = transmissionProbabilityOnContact_noage;

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = relativeTransmissionNoSymptoms;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = riskOfInfectionFromSymptomatic;

    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = recoveredPerInfectedNoSymptoms_noage;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = severePerInfectedSymptoms_noage;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = criticalPerSevere_noage;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = deathsPerCritical_noage;

    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(contact_data_dir, false));
    model.parameters.get<mio::lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.get<mio::lsecir::Seasonality>()     = seasonality;

    auto init = get_initialization(false);
    // Use init as a basis to define appropriate initial values.
    // Compartment values are distributed equally to subcompartments.
    model.populations[0]                   = init[0][0]; //S
    model.populations[LctState::Count - 2] = init[0][6]; //R
    model.populations[LctState::Count - 1] = init[0][7]; //D
    for (size_t i = 1; i < (size_t)InfState::Count - 2; i++) {
        for (size_t subcomp = 0; subcomp < num_subcompartments; subcomp++) {
            model.populations[(i - 1) * num_subcompartments + 1 + subcomp] = init[0][i] / num_subcompartments;
        }
    }

    // Set integrator of fifth order with fixed step size and perform simulation.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max to get a fixed step size.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations = model.calculate_compartments(result);

    if (!save_dir.empty()) {
        std::string filename = save_dir + "lct_nonageresolved_subcomp" + std::to_string(num_subcompartments) + ".h5";
        mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename);
    }

    return mio::success();
}

/** 
* Usage: lct_impact_age_resolution <contact_data_dir> <save_dir> 
*   Both command line arguments are optional but it is beneficial to specify the 
*   <contact_data_dir> as the default is just an educated guess.
*/
int main(int argc, char** argv)
{
    std::string contact_data_dir = "../../data/contacts/";
    std::string save_dir         = "";
    switch (argc) {
    case 3:
        save_dir = argv[2];
        [[fallthrough]];
    case 2:
        contact_data_dir = argv[1];
    }

    ScalarType tmax = 40;

    // Simulation with initial exposed population in age group 2.
    auto result = simulation_with_ageresolution(2, tmax, contact_data_dir, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    // Simulation with initial exposed population in age group 5.
    result = simulation_with_ageresolution(5, tmax, contact_data_dir, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    // Simulation without age resolution.
    result = simulation_without_ageresolution(tmax, contact_data_dir, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
