/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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

#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/io/epi_data.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/parameters.h"
#include "ide_secir/parameters_io.h"
#include "ide_secir/simulation.h"

#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"

#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/filesystem.hpp"
#include <Eigen/src/Core/Random.h>
#include <cstddef>
#include <string>
#include <map>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

// Define parameters for the simulation.
namespace params
{
const size_t num_age_groups        = 6;
const ScalarType age_group_sizes[] = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
const ScalarType total_population  = 83155031.;
// const ScalarType total_confirmed_cases = 341223.;
// const ScalarType deaths = 0.;

// Define transition probabilities per age group.
const ScalarType infectedSymptomsPerInfectedNoSymptoms[] = {0.75, 0.75, 0.8, 0.8, 0.8, 0.8};
const ScalarType severePerInfectedSymptoms[]             = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType criticalPerSevere[]                     = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType deathsPerCritical[]                     = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};

// Define mean stay times per age group.
const ScalarType timeExposed[]            = {4.5, 4.5, 4.5, 4.5, 4.5, 4.5};
const ScalarType timeInfectedNoSymptoms[] = {2.825, 2.825, 2.48, 2.48, 2.48, 2.48};
const ScalarType timeInfectedSymptoms[]   = {7.9895, 7.9895, 7.9734, 7.9139, 7.9139, 7.685};
const ScalarType timeInfectedSevere[]     = {16.855, 16.855, 16.855, 15.61, 13.12, 11.46};
const ScalarType timeInfectedCritical[]   = {17.73, 17.73, 17.064, 17.064, 15.14, 13.66};

// Define number of subcompartments for each compartment per age group.
// These are used as a template argument and have to be constexpr.
constexpr size_t n_subcomps_E[]    = {9, 9, 9, 9, 9, 9};
constexpr size_t n_subcomps_INS[]  = {5, 5, 4, 4, 4, 4};
constexpr size_t n_subcomps_ISy[]  = {16, 16, 16, 15, 15, 13};
constexpr size_t n_subcomps_ISev[] = {8, 8, 8, 7, 6, 5};
constexpr size_t n_subcomps_ICri[] = {8, 8, 8, 8, 7, 6};

// Define lognormal parameters. For each transition, we need shape and scale.
// These are given in this order below. The transition distributions are the same for all age groups.
const ScalarType lognorm_EtINS[]     = {0.32459285, 4.26907484};
const ScalarType lognorm_INStISy[]   = {0.7158751, 0.85135303};
const ScalarType lognorm_INStR[]     = {0.24622068, 7.76114};
const ScalarType lognorm_ISytISev[]  = {0.66258947, 5.29920733};
const ScalarType lognorm_ISytR[]     = {0.24622068, 7.76114};
const ScalarType lognorm_ISevtICri[] = {1.01076765, 0.9};
const ScalarType lognorm_ISevtR[]    = {0.33816427, 17.09411753};
const ScalarType lognorm_ICritD[]    = {0.42819924, 9.76267505};
const ScalarType lognorm_ICritR[]    = {0.33816427, 17.09411753};

// Define epidemiological parameters.
const ScalarType transmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};
const ScalarType relativeTransmissionNoSymptoms     = 1;
const ScalarType riskOfInfectionFromSymptomatic     = 0.3;
const ScalarType seasonality                        = 0.;
const ScalarType scale_confirmed_cases              = 1.;

// Define simulation parameters.
ScalarType t0   = 0.;
ScalarType tmax = 20;
ScalarType dt   = 0.01;

} // namespace params

using namespace mio;

/** 
* @brief Function to transform an age-resolved simulation result into a result without age resolution.
*
* Sums up the values in the age groups to transform the simulation result into a result without age resolution. 
* To provide a clear overview, we use non-age-resolved results for visualizations.
* This implementation is only valid if the simulation is run with equal LctStates for all groups or if the result 
* does not contain any subcompartments (e.g. due to previous accumulation of the subcompartments into compartments).
*   
* @param[in] ageresolved_result TimeSeries with an age-resolved simulation result.
* @returns TimeSeries with the result where the values of the age groups are summed up.
*/
TimeSeries<ScalarType> sum_age_groups(const TimeSeries<ScalarType> ageresolved_result)
{
    using namespace params;
    size_t infstatecount = size_t((ScalarType)ageresolved_result.get_num_elements() / (ScalarType)num_age_groups);
    TimeSeries<ScalarType> nonageresolved_result(infstatecount);

    // For each time point, accumulate the age-resolved result and add the time point
    // to the non-age-resolved result.
    for (Eigen::Index timepoint = 0; timepoint < ageresolved_result.get_num_time_points(); ++timepoint) {
        Eigen::VectorX<ScalarType> result = Eigen::VectorX<ScalarType>::Zero(infstatecount);
        for (size_t infstate = 0; infstate < infstatecount; infstate++) {
            for (size_t group = 0; group < num_age_groups; group++) {
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
IOResult<UncertainContactMatrix<ScalarType>> get_contact_matrix(std::string contact_data_dir,
                                                                bool resolve_by_age = true)
{
    using namespace params;
    const std::string contact_locations[] = {"home", "school_pf_eig", "work", "other"};
    const size_t num_locations            = 4;
    size_t matrix_size                    = num_age_groups;
    if (!resolve_by_age) {
        matrix_size = 1;
    }
    auto contact_matrices = ContactMatrixGroup(num_locations, matrix_size);
    // Load and set baseline contacts for each contact location.
    for (size_t location = 0; location < num_locations; location++) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          read_mobility_plain(contact_data_dir + "baseline_" + contact_locations[location] + ".txt"));
        if (!resolve_by_age) {
            ScalarType average = 0.;
            for (size_t i = 0; i < num_age_groups; i++) {
                for (size_t j = 0; j < num_age_groups; j++) {
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
    return success(UncertainContactMatrix<ScalarType>(contact_matrices));
}

/**
* @brief Function that simulates from time 0 until tmax using an IDE model where we apply a contact scaling after
* two days. The simulation results will be saved in the folder save_dir as .h5 files.    
*
* @param[in] save_dir Directory where simulation results will be saved. Default is an empty string leading to the 
* results not being saved. 
* @returns Any io errors that happen or the simulation results for the compartments.
*/
IOResult<std::pair<Vector, std::vector<std::vector<std::vector<ScalarType>>>>>
simulate_ide(Date start_date, std::string contact_data_dir, std::string reported_data_dir, std::string save_dir = "")
{
    using namespace params;
    using InfTransition = isecir::InfectionTransition;

    // Initialize model.
    // Set total_population_init according to age_group_sizes.
    CustomIndexArray<ScalarType, AgeGroup> total_population_init =
        CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);
    for (size_t group = 0; group < num_age_groups; group++) {
        total_population_init[(AgeGroup)group] = age_group_sizes[group];
    }

    // Set these values to zero since they are set according to RKI data below.
    CustomIndexArray<ScalarType, AgeGroup> deaths_init =
        CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);
    CustomIndexArray<ScalarType, AgeGroup> total_confirmed_cases_init =
        CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), 0.);

    // Initialize with empty series for flows as we will initialize based on RKI data later on.
    isecir::Model model_ide(TimeSeries<ScalarType>((Eigen::Index)isecir::InfectionTransition::Count * num_age_groups),
                            total_population_init, deaths_init, num_age_groups, total_confirmed_cases_init);

    // Set working parameters.
    // Set TransitionDistributions.
    ConstantFunction initialfunc(0);
    StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    std::vector<StateAgeFunctionWrapper> vec_delaydistrib((int)InfTransition::Count, delaydistributioninit);
    // ExposedToInfectedNoSymptoms
    LognormSurvivalFunction survivalExposedToInfectedNoSymptoms(lognorm_EtINS[0], 0, lognorm_EtINS[1]);
    vec_delaydistrib[(int)InfTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
        survivalExposedToInfectedNoSymptoms);
    // InfectedNoSymptomsToInfectedSymptoms
    LognormSurvivalFunction survivalInfectedNoSymptomsToInfectedSymptoms(lognorm_INStISy[0], 0, lognorm_INStISy[1]);
    vec_delaydistrib[(int)InfTransition::InfectedNoSymptomsToInfectedSymptoms].set_state_age_function(
        survivalInfectedNoSymptomsToInfectedSymptoms);
    // InfectedNoSymptomsToRecovered
    LognormSurvivalFunction survivalInfectedNoSymptomsToRecovered(lognorm_INStR[0], 0, lognorm_INStR[1]);
    vec_delaydistrib[(int)InfTransition::InfectedNoSymptomsToRecovered].set_state_age_function(
        survivalInfectedNoSymptomsToRecovered);
    // InfectedSymptomsToInfectedSevere
    LognormSurvivalFunction survivalInfectedSymptomsToInfectedSevere(lognorm_ISytISev[0], 0, lognorm_ISytR[1]);
    vec_delaydistrib[(int)InfTransition::InfectedSymptomsToInfectedSevere].set_state_age_function(
        survivalInfectedSymptomsToInfectedSevere);
    // InfectedSymptomsToRecovered
    LognormSurvivalFunction survivalInfectedSymptomsToRecovered(lognorm_ISytR[0], 0, lognorm_ISytR[1]);
    vec_delaydistrib[(int)InfTransition::InfectedSymptomsToRecovered].set_state_age_function(
        survivalInfectedSymptomsToRecovered);
    // InfectedSevereToInfectedCritical
    LognormSurvivalFunction survivalInfectedSevereToInfectedCritical(lognorm_ISevtICri[0], 0, lognorm_ISevtICri[1]);
    vec_delaydistrib[(int)InfTransition::InfectedSevereToInfectedCritical].set_state_age_function(
        survivalInfectedSevereToInfectedCritical);
    // InfectedSevereToRecovered
    LognormSurvivalFunction survivalInfectedSevereToRecovered(lognorm_ISevtR[0], 0, lognorm_ISevtR[1]);
    vec_delaydistrib[(int)InfTransition::InfectedSevereToRecovered].set_state_age_function(
        survivalInfectedSevereToRecovered);
    // InfectedCriticalToDead
    LognormSurvivalFunction survivalInfectedCriticalToDead(lognorm_ICritD[0], 0, lognorm_ICritR[1]);
    vec_delaydistrib[(int)InfTransition::InfectedCriticalToDead].set_state_age_function(survivalInfectedCriticalToDead);
    // InfectedCriticalToRecovered
    LognormSurvivalFunction survivalInfectedCriticalToRecovered(lognorm_ICritR[0], 0, lognorm_ICritR[1]);
    vec_delaydistrib[(int)InfTransition::InfectedCriticalToRecovered].set_state_age_function(
        survivalInfectedCriticalToRecovered);

    // Set distributions for all age groups.
    for (AgeGroup group = 0; group < (AgeGroup)num_age_groups; group++) {
        model_ide.parameters.get<isecir::TransitionDistributions>()[group] = vec_delaydistrib;
    }

    // Set other parameters.
    for (AgeGroup group = 0; group < (AgeGroup)num_age_groups; group++) {
        std::vector<ScalarType> vec_prob((int)InfTransition::Count, 1.);
        vec_prob[Eigen::Index(InfTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            infectedSymptomsPerInfectedNoSymptoms[(size_t)group];
        vec_prob[Eigen::Index(InfTransition::InfectedNoSymptomsToRecovered)] =
            1 - infectedSymptomsPerInfectedNoSymptoms[(size_t)group];
        vec_prob[Eigen::Index(InfTransition::InfectedSymptomsToInfectedSevere)] =
            severePerInfectedSymptoms[(size_t)group];
        vec_prob[Eigen::Index(InfTransition::InfectedSymptomsToRecovered)] =
            1 - severePerInfectedSymptoms[(size_t)group];
        vec_prob[Eigen::Index(InfTransition::InfectedSevereToInfectedCritical)] = criticalPerSevere[(size_t)group];
        vec_prob[Eigen::Index(InfTransition::InfectedSevereToRecovered)]        = 1 - criticalPerSevere[(size_t)group];
        vec_prob[Eigen::Index(InfTransition::InfectedCriticalToDead)]           = deathsPerCritical[(size_t)group];
        vec_prob[Eigen::Index(InfTransition::InfectedCriticalToRecovered)]      = 1 - deathsPerCritical[(size_t)group];

        model_ide.parameters.get<isecir::TransitionProbabilities>()[group] = vec_prob;
    }

    for (AgeGroup group = 0; group < (AgeGroup)num_age_groups; group++) {
        ConstantFunction constfunc(transmissionProbabilityOnContact[(size_t)group]);
        StateAgeFunctionWrapper StateAgeFunctionWrapperide(constfunc);
        model_ide.parameters.get<isecir::TransmissionProbabilityOnContact>()[group] = StateAgeFunctionWrapperide;

        StateAgeFunctionWrapperide.set_distribution_parameter(relativeTransmissionNoSymptoms);
        model_ide.parameters.get<isecir::RelativeTransmissionNoSymptoms>() = StateAgeFunctionWrapperide;

        StateAgeFunctionWrapperide.set_distribution_parameter(riskOfInfectionFromSymptomatic);
        model_ide.parameters.get<isecir::RiskOfInfectionFromSymptomatic>()[group] = StateAgeFunctionWrapperide;
    }

    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(contact_data_dir, true));
    model_ide.parameters.get<isecir::ContactPatterns>() = contact_matrix;

    model_ide.set_tol_for_support_max(1e-6);

    std::string path_rki = reported_data_dir + "cases_all_age_all_dates.json";
    BOOST_OUTCOME_TRY(std::vector<ConfirmedCasesDataEntry> && rki_data, read_confirmed_cases_data(path_rki));

    // Define vector for scale_confirmed_cases.
    CustomIndexArray<ScalarType, AgeGroup> scale_confirmed_cases_vec =
        CustomIndexArray<ScalarType, AgeGroup>(AgeGroup(num_age_groups), scale_confirmed_cases);

    BOOST_OUTCOME_TRY(isecir::set_initial_flows<ConfirmedCasesDataEntry>(model_ide, dt, rki_data, start_date,
                                                                         scale_confirmed_cases_vec));

    model_ide.check_constraints(dt);

    auto persons_per_state_age = model_ide.get_num_persons_per_state_age(dt);

    // Simulate.
    isecir::Simulation sim(model_ide, dt);
    sim.advance(tmax);

    if (!save_dir.empty()) {
        std::string tmax_string = std::to_string(tmax);
        std::string dt_string   = std::to_string(dt);

        std::string filename_ide = save_dir + "changepoint_ide_" + "_" + tmax_string.substr(0, tmax_string.find(".")) +
                                   "_" + dt_string.substr(0, dt_string.find(".") + 5);

        std::string filename_ide_compartments = filename_ide + "_compartments.h5";
        IOResult<void> save_result_status_c =
            save_result({sim.get_result()}, {0}, num_age_groups, filename_ide_compartments);

        if (!save_result_status_c) {
            return failure(StatusCode::UnknownError, "Error while saving results.");
        }
    }

    // Return init_compartments (i.e. populations at t0) and persons_per_state_age for initialization of other models.

    Vector init_compartments = sim.get_result()[0];

    std::pair<Vector, std::vector<std::vector<std::vector<ScalarType>>>> results =
        std::make_pair(init_compartments, persons_per_state_age);

    // std::vector<TimeSeries<ScalarType>> results = {sim.get_result()};
    return success(results);
}

IOResult<void> simulate_lct(Vector init_compartments, std::string contact_data_dir, std::string save_dir = "")
{
    using namespace params;

    // Initialize age-resolved model.
    using InfState = lsecir::InfectionState;

    using LctState0_4   = LctInfectionState<InfState, 1, n_subcomps_E[0], n_subcomps_INS[0], n_subcomps_ISy[0],
                                            n_subcomps_ISev[0], n_subcomps_ICri[0], 1, 1>;
    using LctState5_14  = LctInfectionState<InfState, 1, n_subcomps_E[1], n_subcomps_INS[1], n_subcomps_ISy[1],
                                            n_subcomps_ISev[1], n_subcomps_ICri[1], 1, 1>;
    using LctState15_34 = LctInfectionState<InfState, 1, n_subcomps_E[2], n_subcomps_INS[2], n_subcomps_ISy[2],
                                            n_subcomps_ISev[2], n_subcomps_ICri[2], 1, 1>;
    using LctState35_59 = LctInfectionState<InfState, 1, n_subcomps_E[3], n_subcomps_INS[3], n_subcomps_ISy[3],
                                            n_subcomps_ISev[3], n_subcomps_ICri[3], 1, 1>;
    using LctState60_79 = LctInfectionState<InfState, 1, n_subcomps_E[4], n_subcomps_INS[4], n_subcomps_ISy[4],
                                            n_subcomps_ISev[4], n_subcomps_ICri[4], 1, 1>;
    using LctState80    = LctInfectionState<InfState, 1, n_subcomps_E[5], n_subcomps_INS[5], n_subcomps_ISy[5],
                                            n_subcomps_ISev[5], n_subcomps_ICri[5], 1, 1>;

    using Model = lsecir::Model<LctState0_4, LctState5_14, LctState15_34, LctState35_59, LctState60_79, LctState80>;

    Model model;

    // Define parameters.
    for (size_t group = 0; group < num_age_groups; group++) {
        model.parameters.get<lsecir::TimeExposed>()[group]            = timeExposed[group];
        model.parameters.get<lsecir::TimeInfectedNoSymptoms>()[group] = timeInfectedNoSymptoms[group];
        model.parameters.get<lsecir::TimeInfectedSymptoms>()[group]   = timeInfectedSymptoms[group];
        model.parameters.get<lsecir::TimeInfectedSevere>()[group]     = timeInfectedSevere[group];
        model.parameters.get<lsecir::TimeInfectedCritical>()[group]   = timeInfectedCritical[group];
        model.parameters.get<lsecir::TransmissionProbabilityOnContact>()[group] =
            transmissionProbabilityOnContact[group];

        model.parameters.get<lsecir::RelativeTransmissionNoSymptoms>()[group] = relativeTransmissionNoSymptoms;
        model.parameters.get<lsecir::RiskOfInfectionFromSymptomatic>()[group] = riskOfInfectionFromSymptomatic;

        model.parameters.get<lsecir::RecoveredPerInfectedNoSymptoms>()[group] =
            1 - infectedSymptomsPerInfectedNoSymptoms[group];
        model.parameters.get<lsecir::SeverePerInfectedSymptoms>()[group] = severePerInfectedSymptoms[group];
        model.parameters.get<lsecir::CriticalPerSevere>()[group]         = criticalPerSevere[group];
        model.parameters.get<lsecir::DeathsPerCritical>()[group]         = deathsPerCritical[group];
    }
    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(contact_data_dir, true));
    model.parameters.get<lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.get<lsecir::Seasonality>()     = seasonality;

    // Use init_compartments as a basis to define appropriate initial values.
    // Compartment values are distributed uniformly to the subcompartments.
    for (size_t group = 0; group < num_age_groups; group++) {
        size_t total_num_subcomps_this_group      = 0;
        size_t total_num_subcomps_previous_groups = 0;
        std::vector<size_t> num_subcompartments;

        switch (group) {
        case 0:
            total_num_subcomps_this_group = LctState0_4::Count;
            num_subcompartments           = {1,
                                             LctState0_4::get_num_subcompartments<InfState::Exposed>(),
                                             LctState0_4::get_num_subcompartments<InfState::InfectedNoSymptoms>(),
                                             LctState0_4::get_num_subcompartments<InfState::InfectedSymptoms>(),
                                             LctState0_4::get_num_subcompartments<InfState::InfectedSevere>(),
                                             LctState0_4::get_num_subcompartments<InfState::InfectedCritical>(),
                                             1,
                                             1};
            break;
        case 1:
            total_num_subcomps_this_group      = LctState5_14::Count;
            total_num_subcomps_previous_groups = LctState0_4::Count;
            num_subcompartments                = {1,
                                                  LctState5_14::get_num_subcompartments<InfState::Exposed>(),
                                                  LctState5_14::get_num_subcompartments<InfState::InfectedNoSymptoms>(),
                                                  LctState5_14::get_num_subcompartments<InfState::InfectedSymptoms>(),
                                                  LctState5_14::get_num_subcompartments<InfState::InfectedSevere>(),
                                                  LctState5_14::get_num_subcompartments<InfState::InfectedCritical>(),
                                                  1,
                                                  1};
            break;
        case 2:
            total_num_subcomps_this_group      = LctState15_34::Count;
            total_num_subcomps_previous_groups = LctState0_4::Count + LctState5_14::Count;
            num_subcompartments                = {1,
                                                  LctState15_34::get_num_subcompartments<InfState::Exposed>(),
                                                  LctState15_34::get_num_subcompartments<InfState::InfectedNoSymptoms>(),
                                                  LctState15_34::get_num_subcompartments<InfState::InfectedSymptoms>(),
                                                  LctState15_34::get_num_subcompartments<InfState::InfectedSevere>(),
                                                  LctState15_34::get_num_subcompartments<InfState::InfectedCritical>(),
                                                  1,
                                                  1};
            break;
        case 3:
            total_num_subcomps_this_group      = LctState35_59::Count;
            total_num_subcomps_previous_groups = LctState0_4::Count + LctState5_14::Count + LctState15_34::Count;
            num_subcompartments                = {1,
                                                  LctState35_59::get_num_subcompartments<InfState::Exposed>(),
                                                  LctState35_59::get_num_subcompartments<InfState::InfectedNoSymptoms>(),
                                                  LctState35_59::get_num_subcompartments<InfState::InfectedSymptoms>(),
                                                  LctState35_59::get_num_subcompartments<InfState::InfectedSevere>(),
                                                  LctState35_59::get_num_subcompartments<InfState::InfectedCritical>(),
                                                  1,
                                                  1};
            break;
        case 4:
            total_num_subcomps_this_group = LctState60_79::Count;
            total_num_subcomps_previous_groups =
                LctState0_4::Count + LctState5_14::Count + LctState15_34::Count + LctState35_59::Count;
            num_subcompartments = {1,
                                   LctState60_79::get_num_subcompartments<InfState::Exposed>(),
                                   LctState60_79::get_num_subcompartments<InfState::InfectedNoSymptoms>(),
                                   LctState60_79::get_num_subcompartments<InfState::InfectedSymptoms>(),
                                   LctState60_79::get_num_subcompartments<InfState::InfectedSevere>(),
                                   LctState60_79::get_num_subcompartments<InfState::InfectedCritical>(),
                                   1,
                                   1};
            break;
        case 5:
            total_num_subcomps_this_group      = LctState80::Count;
            total_num_subcomps_previous_groups = LctState0_4::Count + LctState5_14::Count + LctState15_34::Count +
                                                 LctState35_59::Count + LctState60_79::Count;
            num_subcompartments = {1,
                                   LctState80::get_num_subcompartments<InfState::Exposed>(),
                                   LctState80::get_num_subcompartments<InfState::InfectedNoSymptoms>(),
                                   LctState80::get_num_subcompartments<InfState::InfectedSymptoms>(),
                                   LctState80::get_num_subcompartments<InfState::InfectedSevere>(),
                                   LctState80::get_num_subcompartments<InfState::InfectedCritical>(),
                                   1,
                                   1};
            break;
        }

        model.populations[total_num_subcomps_previous_groups + 0] =
            init_compartments[group * (size_t)isecir::InfectionState::Count + 0]; // Susceptible
        model.populations[total_num_subcomps_previous_groups + total_num_subcomps_this_group - 2] =
            init_compartments[group * (size_t)isecir::InfectionState::Count + 6]; // Recovered
        model.populations[total_num_subcomps_previous_groups + total_num_subcomps_this_group - 1] =
            init_compartments[group * (size_t)isecir::InfectionState::Count + 7]; // Dead
        for (size_t i = (size_t)InfState::Exposed; i < (size_t)InfState::Count - 2; i++) {
            for (size_t subcomp = 0; subcomp < num_subcompartments[i]; subcomp++) {
                model.populations[total_num_subcomps_previous_groups + (i - 1) * num_subcompartments[i] + 1 + subcomp] =
                    init_compartments[group * (size_t)isecir::InfectionState::Count + i] /
                    (ScalarType)num_subcompartments[i];
            }
        }
    }

    // Set integrator of fifth order with fixed step size and perform simulation.
    auto integrator =
        std::make_shared<ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max to get a fixed step size.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    TimeSeries<ScalarType> result = simulate<ScalarType, Model>(0, tmax, dt, model, integrator);
    // Calculate result without division in subcompartments and without division in age groups.
    TimeSeries<ScalarType> populations = sum_age_groups(model.calculate_compartments(result));

    if (!save_dir.empty()) {
        std::string filename              = save_dir + "lct_ageresolved_subcomp" + ".h5";
        IOResult<void> save_result_status = save_result({populations}, {0}, num_age_groups, filename);
    }

    return success();
}

IOResult<void> simulate_ode(Vector init_compartments, std::string contact_data_dir, std::string save_dir = "")
{
    using namespace params;
    // Use ODE FlowModel.
    osecir::Model model_ode(num_age_groups);

    // Set working parameters.
    for (size_t group = 0; group < num_age_groups; group++) {
        model_ode.parameters.get<osecir::TimeExposed<ScalarType>>()[(AgeGroup)group] = timeExposed[group];
        model_ode.parameters.get<osecir::TimeInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group] =
            timeInfectedNoSymptoms[group];
        model_ode.parameters.get<osecir::TimeInfectedSymptoms<ScalarType>>()[(AgeGroup)group] =
            timeInfectedSymptoms[group];
        model_ode.parameters.get<osecir::TimeInfectedSevere<ScalarType>>()[(AgeGroup)group] = timeInfectedSevere[group];
        model_ode.parameters.get<osecir::TimeInfectedCritical<ScalarType>>()[(AgeGroup)group] =
            timeInfectedCritical[group];

        // Set probabilities that determine proportion between compartments.
        model_ode.parameters.get<osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(AgeGroup)group] =
            1 - infectedSymptomsPerInfectedNoSymptoms[group];
        model_ode.parameters.get<osecir::SeverePerInfectedSymptoms<ScalarType>>()[(AgeGroup)group] =
            severePerInfectedSymptoms[group];
        model_ode.parameters.get<osecir::CriticalPerSevere<ScalarType>>()[(AgeGroup)group] = criticalPerSevere[group];
        model_ode.parameters.get<osecir::DeathsPerCritical<ScalarType>>()[(AgeGroup)group] = deathsPerCritical[group];

        // Further model parameters.
        model_ode.parameters.get<osecir::TransmissionProbabilityOnContact<ScalarType>>()[(AgeGroup)group] =
            transmissionProbabilityOnContact[group];
        model_ode.parameters.get<osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[(AgeGroup)group] =
            relativeTransmissionNoSymptoms;
        model_ode.parameters.get<osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(AgeGroup)group] =
            riskOfInfectionFromSymptomatic;
    }
    // Choose TestAndTraceCapacity very large so that riskFromInfectedSymptomatic = RiskOfInfectionFromSymptomatic.
    model_ode.parameters.get<osecir::TestAndTraceCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();
    // Choose ICUCapacity very large so that CriticalPerSevereAdjusted = CriticalPerSevere and deathsPerSevereAdjusted = 0.
    model_ode.parameters.get<osecir::ICUCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();

    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix.
    model_ode.parameters.set<osecir::Seasonality<ScalarType>>(seasonality);

    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(contact_data_dir, true));
    model_ode.parameters.get<osecir::ContactPatterns<ScalarType>>() = contact_matrix;

    // Use  isecir::InfectionState when accessing init_compartments since this is computed using the IDE model.
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::Susceptible}] =
        init_compartments[int(isecir::InfectionState::Susceptible)];
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::Exposed}] =
        init_compartments[int(isecir::InfectionState::Exposed)];
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::InfectedNoSymptoms}] =
        init_compartments[int(isecir::InfectionState::InfectedNoSymptoms)];
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::InfectedSymptoms}] =
        init_compartments[int(isecir::InfectionState::InfectedSymptoms)];
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::InfectedSymptomsConfirmed}] = 0;
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::InfectedSevere}] =
        init_compartments[int(isecir::InfectionState::InfectedSevere)];
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::InfectedCritical}] =
        init_compartments[int(isecir::InfectionState::InfectedCritical)];
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::Recovered}] =
        init_compartments[int(isecir::InfectionState::Recovered)];
    model_ode.populations[{AgeGroup(0), osecir::InfectionState::Dead}] =
        init_compartments[int(isecir::InfectionState::Dead)];

    model_ode.check_constraints();

    // Set integrator and fix step size.
    auto integrator =
        std::make_shared<ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();

    // Simulate.
    std::vector<TimeSeries<ScalarType>> results_ode =
        osecir::simulate_flows<ScalarType>(t0, tmax, dt, model_ode, integrator);

    // Save results.
    if (!save_dir.empty()) {
        std::string tmax_string  = std::to_string(tmax);
        std::string dt_string    = std::to_string(dt);
        std::string filename_ode = save_dir + "ode_" + tmax_string.substr(0, tmax_string.find(".")) + "_" +
                                   dt_string.substr(0, dt_string.find(".") + 5);

        std::string filename_ode_compartments = filename_ode + ".h5";
        IOResult<void> save_result_status_c =
            save_result({results_ode[0]}, {0}, num_age_groups, filename_ode_compartments);
    }

    return success();
}

IOResult<void> simulate_abm(Vector init_compartments, ScalarType tmax, std::string save_dir = "")
{
    unused(init_compartments);
    unused(tmax);
    unused(save_dir);

    return success();
}

int main(int argc, char** argv)
{
    std::string save_dir = "../../simulation_results/";

    // Set result_dir via command line.
    if (argc == 2) {
        save_dir = argv[1];
    }

    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directories(dir);

    Date start_date(2020, 10, 01);

    // Set path to contact data.
    std::string contact_data_dir  = "../../data/Germany/contacts/";
    std::string reported_data_dir = "../../data/Germany/pydata/";

    // Changepoint scenario with halving of contacts after two days.

    auto result_ide = simulate_ide(start_date, contact_data_dir, reported_data_dir, save_dir);

    // Use compartments at time 0 from IDE simulation as initial values for ODE and LCT model to make results comparable.
    auto init_compartments = std::get<0>(result_ide.value());

    auto result_lct = simulate_lct(init_compartments, contact_data_dir, save_dir);

    auto result_ode = simulate_ode(init_compartments, contact_data_dir, save_dir);

    // For use in ABM simulation.
    // The vector has the following structure:
    // - First dimension:  Determines the compartment; we have values for the compartments Expsoed, InfectedNoSymptoms,
    //                     InfectedSymptoms, InfectedSevere and InfectedCritical.
    // - Second dimension: Determines the age group.
    // - Third dimension:  Determines the number of persons per state age. For each element this vector, the
    //                     value yields the number of persons that have a certain state age. The state age can be obtained
    //                     by multiplying the index of the element with the time step size dt.
    std::vector<std::vector<std::vector<ScalarType>>> persons_per_state_age = std::get<1>(result_ide.value());

    // if (!result_ide || !result_lct || !result_ode) {
    //     printf("%s\n", result_ide.error().formatted_message().c_str());
    //     return -1;
    // }

    return 0;
}