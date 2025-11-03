/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler, Julia Bicker
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

#include "abm/infection_state.h"
#include "abm/location_id.h"
#include "abm/location_type.h"
#include "abm/model.h"
#include "abm/parameters.h"
#include "abm/personal_rng.h"
#include "abm/simulation.h"
#include "abm/time.h"
#include "abm/virus_variant.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/io/json_serializer.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/io/epi_data.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
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

#include "config.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/filesystem.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/Random.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <string>
#include <map>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

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
    auto contact_matrices = ContactMatrixGroup<ScalarType>(num_locations, matrix_size);
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
* @brief Function that simulates from time 0 until tmax using an IDE model. The simulation results will be saved in the folder save_dir as .h5 files.    
*
* @param[in] save_dir Directory where simulation results will be saved. Default is an empty string leading to the 
* results not being saved. 
* @returns Any io errors that happen or the simulation results for the compartments.
*/
IOResult<std::pair<Vector, std::vector<std::vector<std::vector<ScalarType>>>>>
simulate_ide(Date start_date, std::string contact_data_dir, std::string reported_data_dir,
             bool exponential_scenario = false, std::string save_dir = "")
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
    // If exponential_scenario is true, we set the transition distributions as exponentially distributed.
    // Else, we use lognormal distributions.

    // Define vector vec_transition_dist_group where transition distributions will be stored.
    ConstantFunction<ScalarType> initialfunc(0);
    StateAgeFunctionWrapper<ScalarType> delaydistributioninit(initialfunc);
    std::vector<StateAgeFunctionWrapper<ScalarType>> vec_transition_dist_group((int)InfTransition::Count,
                                                                               delaydistributioninit);

    if (exponential_scenario) {
        std::cout << "in exp scenario \n";
        // Set TransitionDistributions exponentially distributed.
        for (size_t group = 0; group < num_age_groups; group++) {
            // ExposedToInfectedNoSymptoms
            ExponentialSurvivalFunction<ScalarType> survivalExposedToInfectedNoSymptoms(1. / timeExposed[group]);
            vec_transition_dist_group[(int)InfTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
                survivalExposedToInfectedNoSymptoms);
            // InfectedNoSymptomsToInfectedSymptoms
            ExponentialSurvivalFunction<ScalarType> survivalInfectedNoSymptomsToInfectedSymptoms(
                1. / timeInfectedNoSymptoms[group]);
            vec_transition_dist_group[(int)InfTransition::InfectedNoSymptomsToInfectedSymptoms].set_state_age_function(
                survivalInfectedNoSymptomsToInfectedSymptoms);
            // InfectedNoSymptomsToRecovered
            ExponentialSurvivalFunction<ScalarType> survivalInfectedNoSymptomsToRecovered(
                1. / timeInfectedNoSymptoms[group]);
            vec_transition_dist_group[(int)InfTransition::InfectedNoSymptomsToRecovered].set_state_age_function(
                survivalInfectedNoSymptomsToRecovered);
            // InfectedSymptomsToInfectedSevere
            ExponentialSurvivalFunction<ScalarType> survivalInfectedSymptomsToInfectedSevere(
                1. / timeInfectedSymptoms[group]);
            vec_transition_dist_group[(int)InfTransition::InfectedSymptomsToInfectedSevere].set_state_age_function(
                survivalInfectedSymptomsToInfectedSevere);
            // InfectedSymptomsToRecovered
            ExponentialSurvivalFunction<ScalarType> survivalInfectedSymptomsToRecovered(1. /
                                                                                        timeInfectedSymptoms[group]);
            vec_transition_dist_group[(int)InfTransition::InfectedSymptomsToRecovered].set_state_age_function(
                survivalInfectedSymptomsToRecovered);
            // InfectedSevereToInfectedCritical
            ExponentialSurvivalFunction<ScalarType> survivalInfectedSevereToInfectedCritical(1. /
                                                                                             timeInfectedSevere[group]);
            vec_transition_dist_group[(int)InfTransition::InfectedSevereToInfectedCritical].set_state_age_function(
                survivalInfectedSevereToInfectedCritical);
            // InfectedSevereToRecovered
            ExponentialSurvivalFunction<ScalarType> survivalInfectedSevereToRecovered(1. / timeInfectedSevere[group]);
            vec_transition_dist_group[(int)InfTransition::InfectedSevereToRecovered].set_state_age_function(
                survivalInfectedSevereToRecovered);
            // InfectedCriticalToDead
            ExponentialSurvivalFunction<ScalarType> survivalInfectedCriticalToDead(1. / timeInfectedCritical[group]);
            vec_transition_dist_group[(int)InfTransition::InfectedCriticalToDead].set_state_age_function(
                survivalInfectedCriticalToDead);
            // InfectedCriticalToRecovered
            ExponentialSurvivalFunction<ScalarType> survivalInfectedCriticalToRecovered(1. /
                                                                                        timeInfectedCritical[group]);
            vec_transition_dist_group[(int)InfTransition::InfectedCriticalToRecovered].set_state_age_function(
                survivalInfectedCriticalToRecovered);

            // Set distributions per age group.
            model_ide.parameters.get<isecir::TransitionDistributions>()[(AgeGroup)group] = vec_transition_dist_group;
        }
    }

    else {
        // Set TransitionDistributions lognormally distributed.

        // ExposedToInfectedNoSymptoms
        LognormSurvivalFunction survivalExposedToInfectedNoSymptoms(lognorm_EtINS[0], 0, lognorm_EtINS[1]);
        vec_transition_dist_group[(int)InfTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
            survivalExposedToInfectedNoSymptoms);
        // InfectedNoSymptomsToInfectedSymptoms
        LognormSurvivalFunction survivalInfectedNoSymptomsToInfectedSymptoms(lognorm_INStISy[0], 0, lognorm_INStISy[1]);
        vec_transition_dist_group[(int)InfTransition::InfectedNoSymptomsToInfectedSymptoms].set_state_age_function(
            survivalInfectedNoSymptomsToInfectedSymptoms);
        // InfectedNoSymptomsToRecovered
        LognormSurvivalFunction survivalInfectedNoSymptomsToRecovered(lognorm_INStR[0], 0, lognorm_INStR[1]);
        vec_transition_dist_group[(int)InfTransition::InfectedNoSymptomsToRecovered].set_state_age_function(
            survivalInfectedNoSymptomsToRecovered);
        // InfectedSymptomsToInfectedSevere
        LognormSurvivalFunction survivalInfectedSymptomsToInfectedSevere(lognorm_ISytISev[0], 0, lognorm_ISytR[1]);
        vec_transition_dist_group[(int)InfTransition::InfectedSymptomsToInfectedSevere].set_state_age_function(
            survivalInfectedSymptomsToInfectedSevere);
        // InfectedSymptomsToRecovered
        LognormSurvivalFunction survivalInfectedSymptomsToRecovered(lognorm_ISytR[0], 0, lognorm_ISytR[1]);
        vec_transition_dist_group[(int)InfTransition::InfectedSymptomsToRecovered].set_state_age_function(
            survivalInfectedSymptomsToRecovered);
        // InfectedSevereToInfectedCritical
        LognormSurvivalFunction survivalInfectedSevereToInfectedCritical(lognorm_ISevtICri[0], 0, lognorm_ISevtICri[1]);
        vec_transition_dist_group[(int)InfTransition::InfectedSevereToInfectedCritical].set_state_age_function(
            survivalInfectedSevereToInfectedCritical);
        // InfectedSevereToRecovered
        LognormSurvivalFunction survivalInfectedSevereToRecovered(lognorm_ISevtR[0], 0, lognorm_ISevtR[1]);
        vec_transition_dist_group[(int)InfTransition::InfectedSevereToRecovered].set_state_age_function(
            survivalInfectedSevereToRecovered);
        // InfectedCriticalToDead
        LognormSurvivalFunction survivalInfectedCriticalToDead(lognorm_ICritD[0], 0, lognorm_ICritR[1]);
        vec_transition_dist_group[(int)InfTransition::InfectedCriticalToDead].set_state_age_function(
            survivalInfectedCriticalToDead);
        // InfectedCriticalToRecovered
        LognormSurvivalFunction survivalInfectedCriticalToRecovered(lognorm_ICritR[0], 0, lognorm_ICritR[1]);
        vec_transition_dist_group[(int)InfTransition::InfectedCriticalToRecovered].set_state_age_function(
            survivalInfectedCriticalToRecovered);

        // Set distributions for all age groups since distributions are the same for all groups.
        for (AgeGroup group = 0; group < (AgeGroup)num_age_groups; group++) {
            model_ide.parameters.get<isecir::TransitionDistributions>()[group] = vec_transition_dist_group;
        }
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
        std::string filename_ide = save_dir + "ide" + ".h5";

        // Aggregate age-resolved result to non-age-resolved result.
        TimeSeries<ScalarType> nonageresolved_result = sum_age_groups(sim.get_result());

        // Save non-age-resolved result.
        IOResult<void> save_result_status = save_result({nonageresolved_result}, {0}, 1, filename_ide);
    }

    // Return init_compartments (i.e. populations at t0) and persons_per_state_age for initialization of other models.

    Vector init_compartments = sim.get_result()[0];

    std::pair<Vector, std::vector<std::vector<std::vector<ScalarType>>>> results =
        std::make_pair(init_compartments, persons_per_state_age);

    // std::vector<TimeSeries<ScalarType>> results = {sim.get_result()};
    return success(results);
}

template <bool exponential_scenario>
IOResult<void> simulate_lct(Vector init_compartments, std::string contact_data_dir, std::string save_dir = "")
{
    using namespace params;

    // Initialize age-resolved model.
    using InfState = lsecir::InfectionState;

    using LctState0_4   = LctInfectionState<ScalarType, InfState, 1, n_subcomps_E[0], n_subcomps_INS[0],
                                          n_subcomps_ISy[0], n_subcomps_ISev[0], n_subcomps_ICri[0], 1, 1>;
    using LctState5_14  = LctInfectionState<ScalarType, InfState, 1, n_subcomps_E[1], n_subcomps_INS[1],
                                           n_subcomps_ISy[1], n_subcomps_ISev[1], n_subcomps_ICri[1], 1, 1>;
    using LctState15_34 = LctInfectionState<ScalarType, InfState, 1, n_subcomps_E[2], n_subcomps_INS[2],
                                            n_subcomps_ISy[2], n_subcomps_ISev[2], n_subcomps_ICri[2], 1, 1>;
    using LctState35_59 = LctInfectionState<ScalarType, InfState, 1, n_subcomps_E[3], n_subcomps_INS[3],
                                            n_subcomps_ISy[3], n_subcomps_ISev[3], n_subcomps_ICri[3], 1, 1>;
    using LctState60_79 = LctInfectionState<ScalarType, InfState, 1, n_subcomps_E[4], n_subcomps_INS[4],
                                            n_subcomps_ISy[4], n_subcomps_ISev[4], n_subcomps_ICri[4], 1, 1>;
    using LctState80 = LctInfectionState<ScalarType, InfState, 1, n_subcomps_E[5], n_subcomps_INS[5], n_subcomps_ISy[5],
                                         n_subcomps_ISev[5], n_subcomps_ICri[5], 1, 1>;

    // Define LctState with one subcompartment per compartment for the exponential scenario.
    using LctStateExponential = LctInfectionState<ScalarType, InfState, 1, 1, 1, 1, 1, 1, 1, 1>;

    using Model =
        std::conditional<exponential_scenario,
                         lsecir::Model<ScalarType, LctStateExponential, LctStateExponential, LctStateExponential,
                                       LctStateExponential, LctStateExponential, LctStateExponential>,
                         lsecir::Model<ScalarType, LctState0_4, LctState5_14, LctState15_34, LctState35_59,
                                       LctState60_79, LctState80>>::type;

    Model model;

    // Define parameters.
    for (size_t group = 0; group < num_age_groups; group++) {
        model.parameters.template get<lsecir::TimeExposed<ScalarType>>()[group] = timeExposed[group];
        model.parameters.template get<lsecir::TimeInfectedNoSymptoms<ScalarType>>()[group] =
            timeInfectedNoSymptoms[group];
        model.parameters.template get<lsecir::TimeInfectedSymptoms<ScalarType>>()[group] = timeInfectedSymptoms[group];
        model.parameters.template get<lsecir::TimeInfectedSevere<ScalarType>>()[group]   = timeInfectedSevere[group];
        model.parameters.template get<lsecir::TimeInfectedCritical<ScalarType>>()[group] = timeInfectedCritical[group];
        model.parameters.template get<lsecir::TransmissionProbabilityOnContact<ScalarType>>()[group] =
            transmissionProbabilityOnContact[group];

        model.parameters.template get<lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[group] =
            relativeTransmissionNoSymptoms;
        model.parameters.template get<lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[group] =
            riskOfInfectionFromSymptomatic;

        model.parameters.template get<lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[group] =
            1 - infectedSymptomsPerInfectedNoSymptoms[group];
        model.parameters.template get<lsecir::SeverePerInfectedSymptoms<ScalarType>>()[group] =
            severePerInfectedSymptoms[group];
        model.parameters.template get<lsecir::CriticalPerSevere<ScalarType>>()[group] = criticalPerSevere[group];
        model.parameters.template get<lsecir::DeathsPerCritical<ScalarType>>()[group] = deathsPerCritical[group];
    }
    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(contact_data_dir, true));
    model.parameters.template get<lsecir::ContactPatterns<ScalarType>>() = contact_matrix;
    model.parameters.template get<lsecir::Seasonality<ScalarType>>()     = seasonality;

    // Use init_compartments as a basis to define appropriate initial values.
    // Compartment values are distributed uniformly to the subcompartments.
    for (size_t group = 0; group < num_age_groups; group++) {
        size_t total_num_subcomps_this_group      = 0;
        size_t total_num_subcomps_previous_groups = 0;
        std::vector<size_t> num_subcompartments;

        if (exponential_scenario) {
            total_num_subcomps_this_group      = LctStateExponential::Count;
            total_num_subcomps_previous_groups = group * (size_t)LctStateExponential::Count;
            num_subcompartments                = {1,
                                                  LctStateExponential::get_num_subcompartments<InfState::Exposed>(),
                                                  LctStateExponential::get_num_subcompartments<InfState::InfectedNoSymptoms>(),
                                                  LctStateExponential::get_num_subcompartments<InfState::InfectedSymptoms>(),
                                                  LctStateExponential::get_num_subcompartments<InfState::InfectedSevere>(),
                                                  LctStateExponential::get_num_subcompartments<InfState::InfectedCritical>(),
                                                  1,
                                                  1};
        }

        else {
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
        }

        model.populations[total_num_subcomps_previous_groups + 0] =
            init_compartments[group * (size_t)isecir::InfectionState::Count + 0]; // Susceptible
        model.populations[total_num_subcomps_previous_groups + total_num_subcomps_this_group - 2] =
            init_compartments[group * (size_t)isecir::InfectionState::Count + 6]; // Recovered
        model.populations[total_num_subcomps_previous_groups + total_num_subcomps_this_group - 1] =
            init_compartments[group * (size_t)isecir::InfectionState::Count + 7]; // Dead
        for (size_t i = (size_t)InfState::Exposed; i < (size_t)InfState::Count - 2; i++) {

            // Get the total number of subcomps before the considered InfectionState in this group.
            size_t num_subcomps_previous_infstates =
                std::accumulate(num_subcompartments.begin(), num_subcompartments.begin() + i, 0);

            for (size_t subcomp = 0; subcomp < num_subcompartments[i]; subcomp++) {
                model.populations[total_num_subcomps_previous_groups + num_subcomps_previous_infstates + subcomp] =
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

    if (!save_dir.empty()) {
        std::string filename_lct = save_dir + "lct" + ".h5";

        // Calculate result without division in subcompartments and without division in age groups.
        TimeSeries<ScalarType> nonageresolved_result = sum_age_groups(model.calculate_compartments(result));

        // Save non-age-resolved result.
        IOResult<void> save_result_status = save_result({nonageresolved_result}, {0}, 1, filename_lct);
    }

    return success();
}

IOResult<void> simulate_ode(Vector init_compartments, std::string contact_data_dir, std::string save_dir = "")
{
    using namespace params;
    // Use ODE FlowModel.
    osecir::Model<ScalarType> model_ode(num_age_groups);

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
    for (size_t group = 0; group < num_age_groups; group++) {
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::Susceptible}] =
            init_compartments[group * (size_t)isecir::InfectionState::Count +
                              (size_t)isecir::InfectionState::Susceptible];
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::Exposed}] =
            init_compartments[group * (size_t)isecir::InfectionState::Count + (size_t)isecir::InfectionState::Exposed];
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::InfectedNoSymptoms}] =
            init_compartments[group * (size_t)isecir::InfectionState::Count +
                              (size_t)isecir::InfectionState::InfectedNoSymptoms];
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::InfectedSymptoms}] =
            init_compartments[group * (size_t)isecir::InfectionState::Count +
                              (size_t)isecir::InfectionState::InfectedSymptoms];
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::InfectedSymptomsConfirmed}] = 0;
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::InfectedSevere}] =
            init_compartments[group * (size_t)isecir::InfectionState::Count +
                              (size_t)isecir::InfectionState::InfectedSevere];
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::InfectedCritical}] =
            init_compartments[group * (size_t)isecir::InfectionState::Count +
                              (size_t)isecir::InfectionState::InfectedCritical];
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::Recovered}] =
            init_compartments[group * (size_t)isecir::InfectionState::Count +
                              (size_t)isecir::InfectionState::Recovered];
        model_ode.populations[{AgeGroup(group), osecir::InfectionState::Dead}] =
            init_compartments[group * (size_t)isecir::InfectionState::Count + (size_t)isecir::InfectionState::Dead];
    }

    model_ode.check_constraints();

    // Set integrator and fix step size.
    auto integrator =
        std::make_unique<ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();

    // Simulate.
    std::vector<TimeSeries<ScalarType>> results_ode =
        osecir::simulate_flows<ScalarType>(t0, tmax, dt, model_ode, std::move(integrator));

    // Save results.
    if (!save_dir.empty()) {
        std::string filename_ode = save_dir + "ode" + ".h5";

        // Aggregate age-resolved result to non-age-resolved result.
        TimeSeries<ScalarType> nonageresolved_result = sum_age_groups(results_ode[0]);

        // Save non-age-resolved result.
        IOResult<void> save_result_status = save_result({nonageresolved_result}, {0}, 1, filename_ode);
    }

    return success();
}

/**
 * @brief Aggregate ABM persons' infections states to compartments and save it in time series.
 * @param[in] model Model used for simulation.
 * @param[in] history History containing the time points of the ABM simulation.
 */
mio::TimeSeries<ScalarType> get_abm_compartment_results(
    mio::abm::Model& model,
    mio::History<mio::DataWriterToMemory, ABMLoggers::LogTimePoint, ABMLoggers::LogExposureRate>& history)
{
    // Crate time series. Results are aggregated by infection state and age group.
    mio::TimeSeries<ScalarType> result(static_cast<size_t>(mio::abm::InfectionState::Count) * params::num_age_groups);
    // Time points the ABM simulation made
    auto& tps = std::get<0>(history.get_log());
    for (auto& tp : tps) {
        // Compartments at current time point
        Eigen::VectorXd comps =
            Eigen::VectorXd::Zero(static_cast<size_t>(mio::abm::InfectionState::Count) * params::num_age_groups);
        for (auto& person : model.get_persons()) {
            size_t age = person.get_age().get();
            auto state = person.get_infection_state(tp);
            comps[age * static_cast<size_t>(mio::abm::InfectionState::Count) + static_cast<size_t>(state)] += 1;
        }
        // Time series time points are given in days
        result.add_time_point(tp.days(), comps);
    }
    return result;
}

/**
 * @brief Get flow index from IDE infection transitions.
 * @param[in] old_state Source state of transition.
 * @param[in] new_state Target state of transition.
 */
size_t get_flow_index(mio::abm::InfectionState old_state, mio::abm::InfectionState new_state)
{
    if (old_state == mio::abm::InfectionState::Susceptible) {
        if (new_state == mio::abm::InfectionState::Exposed) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::SusceptibleToExposed);
        }
    }
    else if (old_state == mio::abm::InfectionState::Exposed) {
        if (new_state == mio::abm::InfectionState::InfectedNoSymptoms) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms);
        }
    }
    else if (old_state == mio::abm::InfectionState::InfectedNoSymptoms) {
        if (new_state == mio::abm::InfectionState::InfectedSymptoms) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms);
        }
        else if (new_state == mio::abm::InfectionState::Recovered) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered);
        }
    }
    else if (old_state == mio::abm::InfectionState::InfectedSymptoms) {
        if (new_state == mio::abm::InfectionState::InfectedSevere) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere);
        }
        else if (new_state == mio::abm::InfectionState::Recovered) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered);
        }
    }
    else if (old_state == mio::abm::InfectionState::InfectedSevere) {
        if (new_state == mio::abm::InfectionState::InfectedCritical) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical);
        }
        else if (new_state == mio::abm::InfectionState::Recovered) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::InfectedSevereToRecovered);
        }
    }
    else if (old_state == mio::abm::InfectionState::InfectedCritical) {
        if (new_state == mio::abm::InfectionState::Dead) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::InfectedCriticalToDead);
        }
        else if (new_state == mio::abm::InfectionState::Recovered) {
            return static_cast<size_t>(mio::isecir::InfectionTransition::InfectedCriticalToRecovered);
        }
    }
    return static_cast<size_t>(mio::isecir::InfectionTransition::Count) + 1;
}

/**
 * @brief Compute infection state transitions i.e. number of persons that changed their infection state in the past time step.
 * @param[in] model Model used for simulation.
 * @param[in] history History containing the time points of the ABM simulation.
 */
mio::TimeSeries<ScalarType> get_abm_flows_results(
    mio::abm::Model& model,
    mio::History<mio::DataWriterToMemory, ABMLoggers::LogTimePoint, ABMLoggers::LogExposureRate>& history)
{
    mio::TimeSeries<ScalarType> result(static_cast<size_t>(mio::isecir::InfectionTransition::Count) *
                                       params::num_age_groups);
    auto& tps      = std::get<0>(history.get_log());
    auto tmax      = tps[tps.size() - 1];
    auto t0        = tps[0];
    auto t         = t0;
    auto step_size = mio::abm::hours(1);
    while (t <= tmax) {
        //std::cout << t.days() << "\n";
        // Check if step size is ok i.e. if no person makes more than one transition
        for (auto& person : model.get_persons()) {
            auto current_state = person.get_infection_state(t);
            auto new_state     = person.get_infection_state(t + step_size);
            // Check if person makes transition
            if (current_state != new_state) {
                size_t flow_index = get_flow_index(current_state, new_state);
                if (flow_index == static_cast<size_t>(mio::isecir::InfectionTransition::Count) +
                                      1) { // Person makes more than one transition i.e. the step size has to be adapted
                    auto infection_course = person.get_infection().get_infection_course();
                    int index_current_state;
                    // Find index of current state
                    if (current_state == mio::abm::InfectionState::Susceptible) {
                        index_current_state = -1;
                    }
                    else {
                        auto it = std::find_if(
                            infection_course.begin(), infection_course.end(),
                            [current_state](const std::pair<mio::abm::TimePoint, mio::abm::InfectionState>& a) {
                                return a.second == current_state;
                            });
                        index_current_state = int(std::distance(infection_course.begin(), it));
                    }
                    // Adapt the step size such that it ends one second before the person makes the after next transition
                    auto step_size_new = infection_course[index_current_state + 2].first - t - mio::abm::seconds(1);
                    if (step_size_new >= step_size) {
                        mio::log_error("Error when calculating ABM flow results. Step size has been increased.");
                    }
                    step_size = step_size_new;
                }
            }
        }
        Eigen::VectorXd flows = Eigen::VectorXd::Zero(static_cast<size_t>(mio::isecir::InfectionTransition::Count) *
                                                      params::num_age_groups);
        // Iterate over all persons to compute flows
        for (auto& person : model.get_persons()) {
            size_t age         = person.get_age().get();
            auto current_state = person.get_infection_state(t);
            auto new_state     = person.get_infection_state(t + step_size);
            // Check if person makes transition
            if (current_state != new_state) {
                size_t flow_index = get_flow_index(current_state, new_state);
                // Person makes more than one transition; this should not happen any more because step size was adapted before
                if (flow_index == static_cast<size_t>(mio::isecir::InfectionTransition::Count) + 1) {
                    mio::log_error("No valid infection state transition from {} to {}",
                                   static_cast<size_t>(current_state), static_cast<size_t>(new_state));
                }
                flows[age * static_cast<size_t>(mio::isecir::InfectionTransition::Count) + flow_index] += 1;
            }
        }
        t += step_size;
        // Time in results time series is given in days
        result.add_time_point(t.days(), flows);
        step_size = mio::abm::hours(1);
    }
    return result;
}

/**
 * @brief Calculate the mean transmission rates over time per age group and save it in csv file.
 * @param[in] history History containing the exposure rates per time point.
 * @param[in] save_dir Directory the results are saved in.
 */
void save_transmission_rates(
    mio::History<mio::DataWriterToMemory, ABMLoggers::LogTimePoint, ABMLoggers::LogExposureRate>& history,
    std::string save_dir)
{
    // Calculate cumulative exposure rate over time
    auto& log = std::get<1>(history.get_log());
    std::vector<double> cum_rates(log[0].size());
    for (size_t t = 1; t < log.size(); ++t) {
        for (size_t age = 0; age < params::num_age_groups; age++) {
            cum_rates[age] += log[t][age];
        }
    }

    // Calculate average rate over time
    size_t num_tps = log.size() - 1;
    std::transform(cum_rates.begin(), cum_rates.end(), cum_rates.begin(), [num_tps](double x) {
        return x / num_tps;
    });

    // Save transmission rates in csv file
    std::string filename = save_dir + "transmission_rates.csv";
    auto file            = fopen(filename.c_str(), "w");
    if (file == NULL) {
        mio::log(mio::LogLevel::warn, "Could not open file {}", filename);
    }
    else {
        for (auto rate : cum_rates) {
            fprintf(file, "%.14f,", rate);
        }
        fclose(file);
    }
}

/**
 * @brief Initialize ABM.
 * @param[in] exponential_scenario Specifies whether transition times are exponentially or lognormally distibuted.
 * @param[in] contact_data_dir Directory for contact matrices.
 * @param[in] initial_infection_state_dist Distribution from which initial infection states are drawn.
 * @param[in] rng Rng used for the model and the initialization.
 * @param[in] one_location If true, there is only one location for each location type.
 */
mio::IOResult<mio::abm::Model> initialize_abm(bool exponential_scenario, std::string contact_data_dir,
                                              std::vector<double> initial_infection_state_dist,
                                              mio::RandomNumberGenerator& rng, bool one_location)
{
    using namespace params;
    mio::abm::Model model(num_age_groups);
    model.get_rng() = rng;

    // Set parameters
    if (exponential_scenario) {
        for (size_t group = 0; group < num_age_groups; group++) {
            model.parameters
                .get<mio::abm::TimeExposedToNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] =
                mio::ParameterDistributionExponential(1. / timeExposed[group]);
            model.parameters.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype,
                                                                                mio::AgeGroup(group)}] =
                mio::ParameterDistributionExponential(1. / timeInfectedNoSymptoms[group]);
            model.parameters.get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype,
                                                                                 mio::AgeGroup(group)}] =
                mio::ParameterDistributionExponential(1. / timeInfectedNoSymptoms[group]);
            model.parameters.get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype,
                                                                            mio::AgeGroup(group)}] =
                mio::ParameterDistributionExponential(1. / timeInfectedSymptoms[group]);
            model.parameters.get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype,
                                                                               mio::AgeGroup(group)}] =
                mio::ParameterDistributionExponential(1. / timeInfectedSymptoms[group]);
            model.parameters.get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype,
                                                                            mio::AgeGroup(group)}] =
                mio::ParameterDistributionExponential(1. / timeInfectedSevere[group]);
            model.parameters.get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype,
                                                                             mio::AgeGroup(group)}] =
                mio::ParameterDistributionExponential(1. / timeInfectedSevere[group]);
            model.parameters
                .get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] =
                mio::ParameterDistributionExponential(1. / timeInfectedCritical[group]);
            model.parameters.get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype,
                                                                               mio::AgeGroup(group)}] =
                mio::ParameterDistributionExponential(1. / timeInfectedCritical[group]);
        }
    }
    else { // State transitions are lognormally distributed
        for (size_t group = 0; group < num_age_groups; group++) {
            model.parameters
                .get<mio::abm::TimeExposedToNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] =
                mio::ParameterDistributionLogNormal(lognorm_EtINS[0], lognorm_EtINS[1]);
            model.parameters.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype,
                                                                                mio::AgeGroup(group)}] =
                mio::ParameterDistributionLogNormal(lognorm_INStISy[0], lognorm_INStISy[1]);
            model.parameters.get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype,
                                                                                 mio::AgeGroup(group)}] =
                mio::ParameterDistributionLogNormal(lognorm_INStR[0], lognorm_INStR[1]);
            model.parameters.get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype,
                                                                            mio::AgeGroup(group)}] =
                mio::ParameterDistributionLogNormal(lognorm_ISytISev[0], lognorm_ISytISev[1]);
            model.parameters.get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype,
                                                                               mio::AgeGroup(group)}] =
                mio::ParameterDistributionLogNormal(lognorm_ISytR[0], lognorm_ISytR[1]);
            model.parameters.get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype,
                                                                            mio::AgeGroup(group)}] =
                mio::ParameterDistributionLogNormal(lognorm_ISevtICri[0], lognorm_ISevtICri[1]);
            model.parameters.get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype,
                                                                             mio::AgeGroup(group)}] =
                mio::ParameterDistributionLogNormal(lognorm_ISevtR[0], lognorm_ISevtR[1]);
            model.parameters
                .get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] =
                mio::ParameterDistributionLogNormal(lognorm_ICritD[0], lognorm_ICritD[1]);
            model.parameters.get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype,
                                                                               mio::AgeGroup(group)}] =
                mio::ParameterDistributionLogNormal(lognorm_ICritR[0], lognorm_ICritR[1]);
        }
    }
    for (size_t group = 0; group < num_age_groups; group++) {
        model.parameters
            .get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] =
            infectedSymptomsPerInfectedNoSymptoms[group];
        model.parameters
            .get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] =
            severePerInfectedSymptoms[group];
        model.parameters
            .get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] =
            criticalPerSevere[group];
        model.parameters
            .get<mio::abm::DeathsPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] = 0.;
        model.parameters
            .get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] =
            deathsPerCritical[group];
        // Todo
        model.parameters.get<mio::abm::VirusShedFactor>()[{mio::abm::VirusVariant::Wildtype, mio::AgeGroup(group)}] =
            mio::ParameterDistributionUniform(0., 1.5);
    }

    // We only consider transmission by contacts, therefore aerosol transmission rates are set to 0
    model.parameters.get<mio::abm::AerosolTransmissionRates>() = 0.0;
    // Age group 1 is school-aged
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[mio::AgeGroup(1)] = true;
    // Age groups 2 and 3 are work-aged
    model.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(2)] = true;
    model.parameters.get<mio::abm::AgeGroupGotoWork>()[mio::AgeGroup(3)] = true;

    // Add one hospital and one ICU
    auto hosp_id = model.add_location(mio::abm::LocationType::Hospital);
    auto icu_id  = model.add_location(mio::abm::LocationType::ICU);

    // Map that maps household ids to person ids and number of adults in household
    std::map<mio::abm::LocationId, std::pair<std::vector<mio::abm::PersonId>, size_t>> household_map;

    size_t pop = size_t(total_population);

    // One add one location of each type where all persons are assigned to
    if (one_location) {
        auto home   = model.add_location(mio::abm::LocationType::Home);
        auto work   = model.add_location(mio::abm::LocationType::Work);
        auto school = model.add_location(mio::abm::LocationType::School);
        auto event  = model.add_location(mio::abm::LocationType::SocialEvent);
        auto shop   = model.add_location(mio::abm::LocationType::BasicsShop);
        for (size_t p = 0; p < pop; ++p) {
            auto age = mio::DiscreteDistribution<size_t>::get_instance()(
                model.get_rng(), std::vector<ScalarType>(std::begin(age_group_sizes), std::end(age_group_sizes)));
            // Add person to model
            auto pid     = model.add_person(home, (AgeGroup)age);
            auto& person = model.get_person(pid);
            // Assign locations to person
            person.set_assigned_location(mio::abm::LocationType::Home, home, model.get_id());
            person.set_assigned_location(mio::abm::LocationType::Hospital, hosp_id, model.get_id());
            person.set_assigned_location(mio::abm::LocationType::ICU, icu_id, model.get_id());
            person.set_assigned_location(mio::abm::LocationType::SocialEvent, event, model.get_id());
            person.set_assigned_location(mio::abm::LocationType::BasicsShop, shop, model.get_id());
            if (person.get_age() == mio::AgeGroup(2) || person.get_age() == mio::AgeGroup(3)) { //work-aged
                person.set_assigned_location(mio::abm::LocationType::Work, work, model.get_id());
            }
            if (person.get_age() == mio::AgeGroup(1)) { //school-aged
                person.set_assigned_location(mio::abm::LocationType::School, school, model.get_id());
            }
        }
    }
    else {
        // Add persons
        while (pop > 0) {
            // First create a home location
            auto home = model.add_location(mio::abm::LocationType::Home);
            // Then sample size of home location
            auto home_size =
                mio::DiscreteDistribution<size_t>::get_instance()(
                    model.get_rng(), std::vector<ScalarType>(std::begin(ABMparams::household_size_distribution),
                                                             std::end(ABMparams::household_size_distribution))) +
                1;
            // Check if size is bigger than remaining population
            if (home_size > pop) {
                // If yes, create one household for the remaining population
                home_size = pop;
            }
            // Create persons for the sampled household
            size_t num_adults = 0;
            std::vector<mio::abm::PersonId> person_ids;
            for (size_t p = 0; p < home_size; p++) {
                // Sample person's age group
                auto age = mio::DiscreteDistribution<size_t>::get_instance()(
                    model.get_rng(), std::vector<ScalarType>(std::begin(age_group_sizes), std::end(age_group_sizes)));
                // Add person to model
                auto pid = model.add_person(home, (AgeGroup)age);
                // Add person id to map
                person_ids.push_back(pid);
                // Assign home location to person
                auto& person = model.get_person(pid);
                person.set_assigned_location(mio::abm::LocationType::Home, home, model.get_id());
                // If person has adult age group (AgeGroup >= 2), increase the number of adult counter for this home location
                if (age >= 2) {
                    num_adults++;
                }
            }
            // Add home to household map
            household_map[home] = std::make_pair(person_ids, num_adults);
            // Reduce total population by added persons
            pop -= home_size;
        }

        // Check if all home locations have assigned at least one adult. If there is a home location without adults, another with at least two is searched and one child and one adult are swapped.
        for (auto& [key, value] : household_map) {
            if (value.second == 0) { // value.second is the number of adults
                // Finde home with at least two adults
                auto adult_home_it = std::find_if(household_map.begin(), household_map.end(), [](const auto& pair) {
                    return pair.second.second >= 2;
                });
                if (adult_home_it == household_map.end()) {
                    mio::log_error("No household with more than 1 adult found.");
                    break;
                }
                // Get person ids to change for both households
                auto child_id = value.first[0];
                // Delete one child from children-only-household
                value.first.erase(value.first.begin());
                auto& child   = model.get_person(child_id);
                auto adult_id = mio::abm::PersonId::invalid_ID();
                // Search adult in household with more than one adults
                for (size_t p = 0; p < adult_home_it->second.first.size(); ++p) {
                    auto pid             = adult_home_it->second.first[p];
                    auto& current_person = model.get_person(pid);
                    if (!(current_person.get_age() == mio::AgeGroup(0) ||
                          current_person.get_age() == mio::AgeGroup(1))) {
                        adult_id = pid;
                        // If adult is found, delete adult from its former home location
                        adult_home_it->second.first.erase(adult_home_it->second.first.begin() + p);
                        break;
                    }
                }
                // Get adult from model
                auto& adult = model.get_person(adult_id);
                // Reassign home locations for child and adult
                child.set_assigned_location(mio::abm::LocationType::Home, adult_home_it->first, model.get_id());
                adult.set_assigned_location(mio::abm::LocationType::Home, key, model.get_id());
                // Add child and adult to new homes in the map and increase/decrease num_adults
                household_map[key].first.push_back(adult_id);
                household_map[key].second += 1;
                household_map[adult_home_it->first].first.push_back(child_id);
                household_map[adult_home_it->first].second -= 1;
            }
        }

        // Helper variables for location creation
        size_t curr_work_size     = 0;
        size_t target_work_size   = 0;
        size_t curr_school_size   = 0;
        size_t target_school_size = 0;
        size_t curr_event_size    = 0;
        size_t target_event_size  = 0;
        size_t curr_shop_size     = 0;
        size_t target_shop_size   = 0;

        // LocationIds and number of currently assigned persons per type
        std::vector<std::pair<mio::abm::LocationId, size_t>> works;
        std::vector<std::pair<mio::abm::LocationId, size_t>> schools;
        std::vector<std::pair<mio::abm::LocationId, size_t>> shops;
        std::vector<std::pair<mio::abm::LocationId, size_t>> events;

        // Create locations
        for (auto& person : model.get_persons()) {
            // All agents are assigned hospital and ICU
            person.set_assigned_location(mio::abm::LocationType::Hospital, hosp_id, model.get_id());
            person.set_assigned_location(mio::abm::LocationType::ICU, icu_id, model.get_id());

            // Check if person is work-aged
            if (person.get_age() == mio::AgeGroup(2) || person.get_age() == mio::AgeGroup(3)) {
                // If the current location is full, a new work location has to be created
                if (curr_work_size >= target_work_size) {
                    // Sample work location size
                    size_t size = static_cast<size_t>(mio::NormalDistribution<double>::get_instance()(
                        model.get_rng(), double(ABMparams::workplace_size[0]), double(ABMparams::workplace_size[1])));
                    // Set size to minimum workplace size if necessary
                    target_work_size = std::max({size, ABMparams::min_workplace_size});
                    // Add work location to model and to work location vector
                    auto work_id = model.add_location(mio::abm::LocationType::Work);
                    works.push_back(std::make_pair(work_id, 0));
                    // Set location's capacity to its sampled size
                    auto work = model.get_location(work_id);
                    work.set_capacity(uint32_t(size), uint32_t(1));
                }
                curr_work_size += 1;
            }
            // Check if person is school-aged
            if (person.get_age() == mio::AgeGroup(1)) {
                // If the current school location is full, a new school has to be created
                if (curr_school_size >= target_school_size) {
                    // Sample school size
                    size_t size        = static_cast<size_t>(mio::NormalDistribution<double>::get_instance()(
                        model.get_rng(), double(ABMparams::school_size[0]), double(ABMparams::school_size[1])));
                    target_school_size = std::max({size, ABMparams::min_school_size});
                    // Add school to model and to school location vector
                    auto school_id = model.add_location(mio::abm::LocationType::School);
                    schools.push_back(std::make_pair(school_id, 0));
                    // Set location's capacity to its size
                    auto school = model.get_location(school_id);
                    school.set_capacity(uint32_t(size), uint32_t(1));
                }
                curr_school_size += 1;
            }

            // Event
            if (curr_event_size >= target_event_size) {
                size_t size       = static_cast<size_t>(mio::NormalDistribution<double>::get_instance()(
                    model.get_rng(), double(ABMparams::event_size[0]), double(ABMparams::event_size[1])));
                target_event_size = std::max({size, ABMparams::min_event_size});
                auto event_id     = model.add_location(mio::abm::LocationType::SocialEvent);
                events.push_back(std::make_pair(event_id, 0));
                auto event = model.get_location(event_id);
                event.set_capacity(uint32_t(size), uint32_t(1));
            }
            curr_event_size += 1;
            // Shop
            if (curr_shop_size >= target_shop_size) {
                size_t size      = static_cast<size_t>(mio::NormalDistribution<double>::get_instance()(
                    model.get_rng(), double(ABMparams::shop_size[0]), double(ABMparams::shop_size[1])));
                target_shop_size = std::max({size, ABMparams::min_shop_size});
                auto shop_id     = model.add_location(mio::abm::LocationType::BasicsShop);
                shops.push_back(std::make_pair(shop_id, 0));
                auto shop = model.get_location(shop_id);
                shop.set_capacity(uint32_t(size), uint32_t(1));
            }
            curr_shop_size += 1;
        }
        //Assign locations to persons
        for (auto& person : model.get_persons()) {
            // Check if person is work-aged
            if (person.get_age() == mio::AgeGroup(2) || person.get_age() == mio::AgeGroup(3)) {
                // Sample work uniformly from all work locations
                size_t work_index =
                    mio::UniformIntDistribution<size_t>::get_instance()(model.get_rng(), size_t(0), works.size() - 1);
                person.set_assigned_location(mio::abm::LocationType::Work, works[work_index].first, model.get_id());
                // Increase number of assign persons for sampled work location
                works[work_index].second += 1;
                // Check if work location is full
                if (model.get_location(works[work_index].first).get_capacity().persons <= works[work_index].second) {
                    // If yes, erase location from vector
                    works.erase(works.begin() + work_index);
                }
            }
            // Check if person is school-aged
            if (person.get_age() == mio::AgeGroup(1)) {
                // Sample school uniformly from all school locations
                size_t school_index =
                    mio::UniformIntDistribution<size_t>::get_instance()(model.get_rng(), size_t(0), schools.size() - 1);
                person.set_assigned_location(mio::abm::LocationType::School, schools[school_index].first,
                                             model.get_id());
                // Increase number of assign persons for sampled school locations
                schools[school_index].second += 1;
                // Check if school is full
                if (model.get_location(schools[school_index].first).get_capacity().persons <=
                    schools[school_index].second) {
                    schools.erase(schools.begin() + school_index);
                }
            }
            // Sample event
            size_t event_index =
                mio::UniformIntDistribution<size_t>::get_instance()(model.get_rng(), size_t(0), events.size() - 1);
            person.set_assigned_location(mio::abm::LocationType::SocialEvent, events[event_index].first,
                                         model.get_id());
            events[event_index].second += 1;
            if (model.get_location(events[event_index].first).get_capacity().persons <= events[event_index].second) {
                events.erase(events.begin() + event_index);
            }
            // Sample shop
            size_t shop_index =
                mio::UniformIntDistribution<size_t>::get_instance()(model.get_rng(), size_t(0), shops.size() - 1);
            person.set_assigned_location(mio::abm::LocationType::BasicsShop, shops[shop_index].first, model.get_id());
            shops[shop_index].second += 1;
            if (model.get_location(shops[shop_index].first).get_capacity().persons <= shops[shop_index].second) {
                shops.erase(shops.begin() + shop_index);
            }
        }
    }
    //Set contact rates
    BOOST_OUTCOME_TRY(auto&& home_contacts, read_mobility_plain(contact_data_dir + "baseline_home.txt"));
    BOOST_OUTCOME_TRY(auto&& school_contacts, read_mobility_plain(contact_data_dir + "baseline_school_pf_eig.txt"));
    BOOST_OUTCOME_TRY(auto&& work_contacts, read_mobility_plain(contact_data_dir + "baseline_work.txt"));
    BOOST_OUTCOME_TRY(auto&& other_contacts, read_mobility_plain(contact_data_dir + "baseline_other.txt"));

    for (auto& loc : model.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            for (size_t from = 0; from < num_age_groups; from++) {
                for (size_t to = 0; to < num_age_groups; to++) {
                    loc.get_infection_parameters()
                        .get<mio::abm::ContactRates>()[{mio::AgeGroup(from), mio::AgeGroup(to)}] =
                        home_contacts(from, to);
                    ;
                }
            }
            break;
        case mio::abm::LocationType::School:
            for (size_t from = 0; from < num_age_groups; from++) {
                for (size_t to = 0; to < num_age_groups; to++) {
                    loc.get_infection_parameters()
                        .get<mio::abm::ContactRates>()[{mio::AgeGroup(from), mio::AgeGroup(to)}] =
                        school_contacts(from, to);
                }
            }
            break;
        case mio::abm::LocationType::Work:
            for (size_t from = 0; from < num_age_groups; from++) {
                for (size_t to = 0; to < num_age_groups; to++) {
                    loc.get_infection_parameters()
                        .get<mio::abm::ContactRates>()[{mio::AgeGroup(from), mio::AgeGroup(to)}] =
                        work_contacts(from, to);
                }
            }
            break;
        case mio::abm::LocationType::SocialEvent:
            for (size_t from = 0; from < num_age_groups; from++) {
                for (size_t to = 0; to < num_age_groups; to++) {
                    loc.get_infection_parameters()
                        .get<mio::abm::ContactRates>()[{mio::AgeGroup(from), mio::AgeGroup(to)}] =
                        other_contacts(from, to);
                }
            }
            break;
        case mio::abm::LocationType::BasicsShop:
            for (size_t from = 0; from < num_age_groups; from++) {
                for (size_t to = 0; to < num_age_groups; to++) {
                    loc.get_infection_parameters()
                        .get<mio::abm::ContactRates>()[{mio::AgeGroup(from), mio::AgeGroup(to)}] =
                        other_contacts(from, to);
                }
            }
            break;
        default:
            for (size_t from = 0; from < num_age_groups; from++) {
                for (size_t to = 0; to < num_age_groups; to++) {
                    loc.get_infection_parameters()
                        .get<mio::abm::ContactRates>()[{mio::AgeGroup(from), mio::AgeGroup(to)}] = 1.;
                }
            }
            break;
        }
    }

    // Initialize infection states
    for (auto& person : model.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(model.get_rng(), initial_infection_state_dist));
        auto p_rng = mio::abm::PersonalRandomNumberGenerator(person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(p_rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, mio::abm::TimePoint(int(t0)),
                                                         infection_state));
        }
    }

    return mio::success(model);
}

/**
 * @brief Run one ABM simulation.
 * @param[in] exponential_scenario Specifies whether transition times are exponentially or lognormally distibuted.
 * @param[in] tmax End time of simulation.
 * @param[in] contact_data_dir Directory for contact matrices.
 * @param[in] initial_infection_state_dist Distribution from which initial infection states are drawn.
 * @param[in] save_results If true the compartment time series, flow time series, the model and the average exposure rates over time are saved after the simulation.
 * @param[in] save_dir Directory the results are saved in and/or from which the model is read in.
 * @param[in] read_model Spacifies whether the model is read in or is initialized from scratch with initialize_abm. If the model is read in, the simulation is done using the given rng.
 * @param[in] rng RNG used for model simulation.
 * @param[in] start_time Start time of the simulation in days.
 * @param[in] one_location If true, all agents there is only one location for each type.
 */
IOResult<mio::TimeSeries<double>> simulate_abm(bool exponential_scenario, ScalarType tmax, std::string contact_data_dir,
                                               std::vector<double> initial_infection_state_dist, bool save_results,
                                               std::string save_dir, bool read_model, mio::RandomNumberGenerator& rng,
                                               double start_time, bool one_location)
{
    mio::set_log_level(mio::LogLevel::warn);
    mio::abm::TimePoint start_tp = mio::abm::TimePoint(mio::abm::days(int(start_time)).seconds());
    using namespace params;
    mio::abm::Model model(num_age_groups);
    if (read_model) {
        std::cout << "Reading in and initializing model...\n";
        BOOST_OUTCOME_TRY(auto&& model_for_initialization,
                          mio::read_json(save_dir + "model.json", mio::Tag<mio::abm::Model>()));
        // Set model rng
        model.get_rng() = rng;
        // Set parameters
        model.parameters = model_for_initialization.parameters;
        // Copy locations
        for (auto&& location : model_for_initialization.get_locations()) {
            if (location.get_type() != mio::abm::LocationType::Cemetery) {
                model.add_location(location.get_type());
            }
        }
        // Add persons with assigned locations and infection states
        for (auto&& person : model_for_initialization.get_persons()) {
            auto home_id         = person.get_assigned_location(mio::abm::LocationType::Home);
            auto infection_state = person.get_infection_state(start_tp);
            auto& new_person     = model.get_person(model.add_person(home_id, person.get_age()));
            // Set assigned locations
            new_person.set_assigned_location(mio::abm::LocationType::Home, home_id, model.get_id());
            new_person.set_assigned_location(mio::abm::LocationType::Hospital,
                                             person.get_assigned_location(mio::abm::LocationType::Hospital),
                                             model.get_id());
            new_person.set_assigned_location(mio::abm::LocationType::ICU,
                                             person.get_assigned_location(mio::abm::LocationType::ICU), model.get_id());
            new_person.set_assigned_location(mio::abm::LocationType::SocialEvent,
                                             person.get_assigned_location(mio::abm::LocationType::SocialEvent),
                                             model.get_id());
            new_person.set_assigned_location(mio::abm::LocationType::BasicsShop,
                                             person.get_assigned_location(mio::abm::LocationType::BasicsShop),
                                             model.get_id());
            if (person.get_age() == mio::AgeGroup(2) || person.get_age() == mio::AgeGroup(3)) {
                new_person.set_assigned_location(mio::abm::LocationType::Work,
                                                 person.get_assigned_location(mio::abm::LocationType::Work),
                                                 model.get_id());
            }
            else if (person.get_age() == mio::AgeGroup(1)) {
                new_person.set_assigned_location(mio::abm::LocationType::School,
                                                 person.get_assigned_location(mio::abm::LocationType::School),
                                                 model.get_id());
            }

            // Set infection state
            if (infection_state != mio::abm::InfectionState::Susceptible) {
                auto p_rng = mio::abm::PersonalRandomNumberGenerator(new_person);
                new_person.add_new_infection(mio::abm::Infection(p_rng, mio::abm::VirusVariant::Wildtype,
                                                                 new_person.get_age(), model.parameters, start_tp,
                                                                 infection_state));
            }
        }
        std::cout << "Reading in and initializing model finished!\n";
    }
    else {
        std::cout << "Initializing model...\n";
        BOOST_OUTCOME_TRY(model, initialize_abm(exponential_scenario, contact_data_dir, initial_infection_state_dist,
                                                rng, one_location));
        std::cout << "Initializing model finished!\n";
    }

    mio::History<mio::DataWriterToMemory, ABMLoggers::LogTimePoint, ABMLoggers::LogExposureRate> history;
    auto sim = mio::abm::Simulation(start_tp, std::move(model));
    std::cout << "Advancing the simulation...\n";
    sim.advance(start_tp + mio::abm::days(int(tmax)), history);
    std::cout << "Simulation finished!\n";

    auto compartment_result = get_abm_compartment_results(sim.get_model(), history);
    if (save_results) {
        // Save compartments
        std::cout << "Exporting compartment results...\n";
        BOOST_OUTCOME_TRY(compartment_result.export_csv(save_dir + "comps.csv"))
        // Save flows
        std::cout << "Computing result flows...\n";
        auto flow_result = get_abm_flows_results(sim.get_model(), history);
        std::cout << "Exporting flows...\n";
        BOOST_OUTCOME_TRY(flow_result.export_csv(save_dir + "flows.csv"))
        // Save model state
        std::cout << "Saving model in json file...\n";
        BOOST_OUTCOME_TRY(mio::write_json(save_dir + "model.json", sim.get_model()));
        std::cout << "Saving transmission rates...\n";
        save_transmission_rates(history, save_dir);
        std::cout << "Saving results finished!\n";
    }

    return success(compartment_result);
}

/**
 * @brief Perform ABM ensemble run.
 * @param[in] num_runs Number of runs performed for the ensemble.
 * @param[in] exponential_scenario Specifies whether the state transition times are exponentially or lognormally distributed.
 * @param[in] tmax End time of simulation.
 * @param[in] contact_data_dir Directory for contact matrices.
 * @param[in] initial_infection_state_dist Initial infection state distribution. Only used if model is not read in.
 * @param[in] save_dir Output directory.
 * @param[in] read_model Specifies whether the model should be read in for each simulation.
 * @param[in] save_single_results Specifies whether all individual simulation results should be saved.
 * @param[in] start_time Start time of the simulation in days.
 * @param[in] one_location If true, there is only one location of each location type.
 */
mio::IOResult<void> abm_ensemble_run(size_t num_runs, bool exponential_scenario, ScalarType tmax,
                                     std::string contact_data_dir, std::vector<double> initial_infection_state_dist,
                                     std::string save_dir, bool read_model, bool save_single_results, double start_time,
                                     bool one_location)
{
    std::vector<std::vector<mio::TimeSeries<double>>> ensemble_result(
        num_runs, std::vector<mio::TimeSeries<ScalarType>>(
                      1, mio::TimeSeries<ScalarType>(int(mio::abm::InfectionState::Count) * params::num_age_groups)));
    for (size_t run = 0; run < num_runs; ++run) {
        std::cout << "Running run " << std::to_string(run) << " of " << std::to_string(num_runs) << std::endl;
        // Generate rng for current run
        mio::RandomNumberGenerator rng;
        //Run ABM simulation
        BOOST_OUTCOME_TRY(ensemble_result[run][0],
                          simulate_abm(exponential_scenario, tmax, contact_data_dir, initial_infection_state_dist,
                                       save_single_results, save_dir, read_model, rng, start_time, one_location));
    }
    // Calculate ensemble percentiles
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.05)[0].export_csv(save_dir + "ABM_p05.csv"));
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.25)[0].export_csv(save_dir + "ABM_p25.csv"));
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.50)[0].export_csv(save_dir + "ABM_p50.csv"));
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.75)[0].export_csv(save_dir + "ABM_p75.csv"));
    BOOST_OUTCOME_TRY(mio::ensemble_percentile(ensemble_result, 0.95)[0].export_csv(save_dir + "ABM_p95.csv"));
    return mio::success();
}

int main()
{
    constexpr bool exponential_scenario = true;
    const std::vector<double> infection_distribution{0.99, 0.005, 0.005, 0.0, 0.0, 0.0, 0.0, 0.0};

    std::string save_dir = "";
    if (exponential_scenario) {
        save_dir = "../../simulation_results/compare_abm_ide_lct_ode/exponential/";
    }
    else {
        save_dir = "../../simulation_results/compare_abm_ide_lct_ode/different_dists/";
    }

    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directories(dir);

    Date start_date(2020, 10, 01);

    // Set path to contact data.
    std::string contact_data_dir  = "../../data/Germany/contacts/";
    std::string reported_data_dir = "../../data/Germany/pydata/";

    // Simulate one ABM run to get initialization for following simulations
    std::vector<uint32_t> seeds = {518254298, 179139073, 1937166324, 3882038653, 1776323092, 1261445412};
    std::string set_name        = "Seed1";
    save_dir += set_name + "/";
    boost::filesystem::path abm_dir(save_dir);
    boost::filesystem::create_directories(abm_dir);
    mio::RandomNumberGenerator rng;
    rng.seed(seeds);
    double init_tmax  = 14;
    bool one_location = false;
    // Simulation for initialization is 14 days
    auto init_result = simulate_abm(exponential_scenario, init_tmax, contact_data_dir, infection_distribution, true,
                                    save_dir, false, rng, params::t0, one_location);
    // Simulate ABM ensemble run.
    size_t num_runs = 100;
    auto result_abm =
        abm_ensemble_run(num_runs, exponential_scenario, params::t0 + init_tmax + params::tmax, contact_data_dir,
                         infection_distribution, save_dir, true, false, params::t0 + init_tmax, one_location);

    // Simulate IDE.
    // auto result_ide = simulate_ide(start_date, contact_data_dir, reported_data_dir, exponential_scenario, save_dir);

    // // Use compartments at time 0 from IDE simulation as initial values for ODE and LCT model to make results comparable.
    // auto init_compartments = std::get<0>(result_ide.value());

    // // Simulate LCT.
    // auto result_lct = simulate_lct<exponential_scenario>(init_compartments, contact_data_dir, save_dir);

    // // Simulate ODE.
    // auto result_ode = simulate_ode(init_compartments, contact_data_dir, save_dir);

    // For use in ABM simulation.
    // The vector has the following structure:
    // - First dimension:  Determines the compartment; we have values for the compartments Expsoed, InfectedNoSymptoms,
    //                     InfectedSymptoms, InfectedSevere and InfectedCritical.
    // - Second dimension: Determines the age group.
    // - Third dimension:  Determines the number of persons per state age. For each element this vector, the
    //                     value yields the number of persons that have a certain state age. The state age can be obtained
    //                     by multiplying the index of the element with the time step size dt.
    //std::vector<std::vector<std::vector<ScalarType>>> persons_per_state_age = std::get<1>(result_ide.value());

    // if (!result_ide || !result_lct || !result_ode) {
    //     printf("%s\n", result_ide.error().formatted_message().c_str());
    //     return -1;
    // }

    return 0;
}
