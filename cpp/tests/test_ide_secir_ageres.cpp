/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Hannah Tritzschak, Anna Wendler
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
#include "load_test_data.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <gtest/gtest.h>
#include <vector>

TEST(TestIdeAgeres, compareWithPreviousRun)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    size_t num_agegroups = 3;
    ScalarType tmax      = 5.0;
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> N =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 5000.);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 6.);
    ScalarType dt = 1.;

    int num_transitions  = (int)mio::isecir::InfectionTransition::Count;
    int num_compartments = (int)mio::isecir::InfectionState::Count;

    // Create TimeSeries with num_transitions * num_agegroups elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_transitions * num_agegroups);

    // Define transitions that will be used for initialization.
    Vec vec_init = Vec::Constant(num_transitions * num_agegroups, 1.);
    // First AgeGroup.
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 20.0;
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 3.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 10.0;
    // Second AgeGroup.
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 20.0;
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
    vec_init[num_transitions + (int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 10.0;
    // Third AgeGroup.
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 25.0;
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
    vec_init[2 * num_transitions + (int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
    // Add initial time point to TimeSeries.
    init.add_time_point(-10., vec_init);
    // Add further time points until t0.
    while (init.get_last_time() < -dt / 2.) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, deaths, num_agegroups);

    // Set working parameters.

    // First AgeGroup for TransitionDistributions.
    mio::SmootherCosine<ScalarType> smoothcos1(2.0);
    mio::StateAgeFunctionWrapper<ScalarType> delaydistribution1(smoothcos1);
    std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_delaydistrib1(num_transitions, delaydistribution1);
    // TransitionDistribution is not used for SusceptibleToExposed. Therefore, the parameter can be set to any value.
    vec_delaydistrib1[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_distribution_parameter(-1.);
    vec_delaydistrib1[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_distribution_parameter(3.0);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(0)] = vec_delaydistrib1;

    // Second AgeGroup for TransitionDistributions.
    mio::SmootherCosine<ScalarType> smoothcos2(3.0);
    mio::StateAgeFunctionWrapper<ScalarType> delaydistribution2(smoothcos2);
    std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_delaydistrib2(num_transitions, delaydistribution2);
    // TransitionDistribution is not used for SusceptibleToExposed. Therefore, the parameter can be set to any value.
    vec_delaydistrib2[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_distribution_parameter(-1.);
    vec_delaydistrib2[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_distribution_parameter(2.0);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(1)] = vec_delaydistrib2;

    // Third AgeGroup for TransitionDistributions.
    mio::SmootherCosine<ScalarType> smoothcos3(2.5);
    mio::StateAgeFunctionWrapper<ScalarType> delaydistribution3(smoothcos3);
    std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_delaydistrib3(num_transitions, delaydistribution3);
    // TransitionDistribution is not used for SusceptibleToExposed. Therefore, the parameter can be set to any value.
    vec_delaydistrib1[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_distribution_parameter(-1.);
    vec_delaydistrib1[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_distribution_parameter(4.0);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(2)] = vec_delaydistrib3;

    //TransitionProbabilities.
    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        std::vector<ScalarType> vec_prob(num_transitions, 0.5);
        // The following probabilities must be 1, as there is no other way to go.
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        model.parameters.get<mio::isecir::TransitionProbabilities>()[group]                   = vec_prob;
    }

    // Contact matrix.
    mio::ContactMatrixGroup<ScalarType> contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, num_agegroups);
    contact_matrix[0] =
        mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(num_agegroups, num_agegroups, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ExponentialSurvivalFunction exponential(0.5);
    mio::StateAgeFunctionWrapper<ScalarType> prob(exponential);
    for (auto group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        model.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[group] = prob;
        model.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[group]   = prob;
        model.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[group]   = prob;
    }
    model.parameters.set<mio::isecir::Seasonality>(0.1);
    // Start the simulation on the 40th day of a year (i.e. in February).
    model.parameters.set<mio::isecir::StartDay>(40);

    model.check_constraints(dt);

    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);

    // Compare compartments at last time point with results from a previous run that are given here.
    mio::TimeSeries<ScalarType> compartments = sim.get_result();
    Eigen::VectorX<ScalarType> compare_compartments(num_compartments * num_agegroups);
    compare_compartments << 469.6949422507, 18.2547785569, 26.0031050672, 7.9621653957, 3.7350284051, 1.7375751388,
        4461.7342012617, 10.8782039239, 469.6949422507, 35.8652960540, 24.7473405359, 15.4322222059, 6.7061736912,
        2.9002347700, 4434.3772098068, 10.2765806855, 587.1186778134, 33.9201025102, 31.8752920696, 14.6856309847,
        6.6453275434, 2.9575629755, 4311.2552183628, 11.5421877404;

    ASSERT_EQ(compare_compartments.size(), static_cast<size_t>(compartments.get_last_value().size()));

    for (int j = 0; j < compare_compartments.size(); j++) {
        ASSERT_NEAR(compartments.get_last_value()[j], compare_compartments[j], 1e-7);
    }

    // Compare transitions at last time point with results from a previous run that are given here.

    mio::TimeSeries<ScalarType> transitions = sim.get_transitions();
    Eigen::VectorX<ScalarType> compare_transitions(num_transitions * num_agegroups);
    compare_transitions << 36.5095571138, 35.2210349942, 15.9243307913, 16.7851751402, 7.4700568103, 7.4700568103,
        3.4751502776, 3.4751502776, 1.5959804568, 1.5959804568, 36.5095571138, 33.5703502804, 15.9243307913,
        14.9401136206, 6.9503005552, 6.9503005552, 3.0041079341, 3.0041079341, 1.3063243113, 1.3063243113,
        45.6369463923, 43.0480501661, 19.8862347266, 19.8862347266, 9.0259862757, 9.0259862757, 4.0260421953,
        4.0260421953, 1.7699583766, 1.7699583766;

    ASSERT_EQ(compare_transitions.size(), static_cast<size_t>(transitions.get_last_value().size()));

    for (int j = 0; j < compare_transitions.size(); j++) {
        ASSERT_NEAR(transitions.get_last_value()[j], compare_transitions[j], 1e-7);
    }
}

// Check results of a simulation with two age groups with an example calculated by hand,
// see https://doi.org/10.1016/j.amc.2025.129636 for the used formulas.
// The idea of this test is to mimic the checkSimulationFunctions test for the nonageresolved IDE-SECIR model.
// For this, we set up two models with two age groups, respectively. All parameters are set as in the non-ageresolved
// test. In the first model, we only have individuals in the first age group where the population size as well as the
// initial values for the deaths and transitions are chosen as in the non-ageresolved test. Hence, we can use the
// values from the other test to check for correctness. In the second model, we only have individuals in the second
// age group and the population size as well as the initial values for the deaths and transitions are scaled by some
// factor. This way, the expected results for this model are determined by the results from the other test and the
// scaling factor. With this, we can check that the correct indices are used for each age group.
TEST(IdeSecirAgeres, checkSimulationFunctions)
{
    using Vec            = mio::TimeSeries<ScalarType>::Vector;
    size_t num_agegroups = 2;
    ScalarType tmax      = 0.5;
    ScalarType dt        = 0.5;

    // Define baseline values for the population size, the number of deaths and values for the transitions
    // from SusceptibleToExposed and InfectedNoSymptomsToInfectedSymptoms for the considered age group.
    // These values correspond to the values of the non-ageresolved checkSimulationFunctions test.
    ScalarType population_baseline                           = 10000.;
    ScalarType deaths_baseline                               = 10.;
    ScalarType susceptibleToExposed_baseline                 = 1.;
    ScalarType infectedNoSymptomsToInfectedSymptoms_baseline = 8.;

    // Scaling factor that determines by which factor the baseline values for the initialization will be scaled for
    // the second age group.
    ScalarType baseline_scaling = 2.;

    /*****************************************************************************************************************/

    // Model with whole population in first age group with index 0. Use baseline values.
    int considered_group                                      = 0;
    std::vector<ScalarType> N_vec_group0                      = {population_baseline, 0};
    std::vector<ScalarType> deaths_vec_group0                 = {deaths_baseline, 0};
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> N_group0 = mio::CustomIndexArray<ScalarType, mio::AgeGroup>(
        mio::AgeGroup(num_agegroups), N_vec_group0.begin(), N_vec_group0.end());
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths_group0 = mio::CustomIndexArray<ScalarType, mio::AgeGroup>(
        mio::AgeGroup(num_agegroups), deaths_vec_group0.begin(), deaths_vec_group0.end());

    // Create TimeSeries with num_transitions elements where transitions needed for simulation will be stored.
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> init_group0(num_transitions * num_agegroups);

    // Add time points for initialization for transitions and death.
    Vec vec_init_group0 = Vec::Constant(num_transitions * num_agegroups, 0.);
    // Set values for vec_init for group 0 as in nonageresolved test.
    vec_init_group0[considered_group * (int)mio::isecir::InfectionTransition::Count +
                    (int)mio::isecir::InfectionTransition::SusceptibleToExposed] = susceptibleToExposed_baseline;
    vec_init_group0[considered_group * (int)mio::isecir::InfectionTransition::Count +
                    (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
        infectedNoSymptomsToInfectedSymptoms_baseline;

    // Add initial time point to TimeSeries.
    init_group0.add_time_point(-0.5, vec_init_group0);
    while (init_group0.get_last_time() < 0) {
        init_group0.add_time_point(init_group0.get_last_time() + dt, vec_init_group0);
    }

    // Initialize model.
    mio::isecir::Model model_group0(std::move(init_group0), N_group0, deaths_group0, num_agegroups);

    // Set working parameters.
    for (size_t group = 0; group < num_agegroups; group++) {
        // In this example, SmootherCosine with parameter 1 (and thus with a maximum support of 1)
        // is used for all TransitionDistribution%s for all groups.
        mio::SmootherCosine<ScalarType> smoothcos(1.0);
        mio::StateAgeFunctionWrapper<ScalarType> delaydistribution(smoothcos);
        std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_delaydistrib(num_transitions, delaydistribution);
        model_group0.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(group)] = vec_delaydistrib;

        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]            = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]     = 1;
        model_group0.parameters.get<mio::isecir::TransitionProbabilities>()[mio::AgeGroup(group)] = vec_prob;

        mio::ContactMatrixGroup<ScalarType> contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, num_agegroups);
        contact_matrix[0] =
            mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(num_agegroups, num_agegroups, 4.));
        model_group0.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        mio::SmootherCosine<ScalarType> smoothcos_prob(1.0);
        mio::StateAgeFunctionWrapper<ScalarType> prob(smoothcos_prob);
        model_group0.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[mio::AgeGroup(group)] = prob;
        model_group0.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[mio::AgeGroup(group)]   = prob;
        model_group0.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[mio::AgeGroup(group)]   = prob;
    }
    model_group0.parameters.set<mio::isecir::Seasonality>(0.);
    model_group0.parameters.set<mio::isecir::StartDay>(0);

    // Carry out simulation.
    mio::isecir::Simulation sim_group0(model_group0, dt);
    sim_group0.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated_group0    = sim_group0.get_result();
    mio::TimeSeries<ScalarType> transitions_simulated_group0 = sim_group0.get_transitions();

    /*****************************************************************************************************************/

    // Model with whole population in second age group with index 1. Use baseline values scaled by baseline_scaling.
    considered_group                                          = 1;
    std::vector<ScalarType> N_vec_group1                      = {0., baseline_scaling * population_baseline};
    std::vector<ScalarType> deaths_vec_group1                 = {0., baseline_scaling * deaths_baseline};
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> N_group1 = mio::CustomIndexArray<ScalarType, mio::AgeGroup>(
        mio::AgeGroup(num_agegroups), N_vec_group1.begin(), N_vec_group1.end());
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths_group1 = mio::CustomIndexArray<ScalarType, mio::AgeGroup>(
        mio::AgeGroup(num_agegroups), deaths_vec_group1.begin(), deaths_vec_group1.end());

    // Create TimeSeries with num_transitions elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init_group1(num_transitions * num_agegroups);

    // Add time points for initialization for transitions and death.
    Vec vec_init_group1 = Vec::Constant(num_transitions * num_agegroups, 0.);
    // Set values for vec_init for group 0 as in nonageresolved test.
    vec_init_group1[considered_group * (int)mio::isecir::InfectionTransition::Count +
                    (int)mio::isecir::InfectionTransition::SusceptibleToExposed] =
        baseline_scaling * susceptibleToExposed_baseline;
    vec_init_group1[considered_group * (int)mio::isecir::InfectionTransition::Count +
                    (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
        baseline_scaling * infectedNoSymptomsToInfectedSymptoms_baseline;

    // Add initial time point to TimeSeries.
    init_group1.add_time_point(-0.5, vec_init_group1);
    while (init_group1.get_last_time() < 0) {
        init_group1.add_time_point(init_group1.get_last_time() + dt, vec_init_group1);
    }

    // Initialize model.
    mio::isecir::Model model_group1(std::move(init_group1), N_group1, deaths_group1, num_agegroups);

    // Set working parameters.
    for (size_t group = 0; group < num_agegroups; group++) {
        // In this example, SmootherCosine with parameter 1 (and thus with a maximum support of 1)
        // is used for all TransitionDistribution%s for all groups.
        mio::SmootherCosine<ScalarType> smoothcos(1.0);
        mio::StateAgeFunctionWrapper<ScalarType> delaydistribution(smoothcos);
        std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_delaydistrib(num_transitions, delaydistribution);
        model_group1.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(group)] = vec_delaydistrib;

        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]            = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]     = 1;
        model_group1.parameters.get<mio::isecir::TransitionProbabilities>()[mio::AgeGroup(group)] = vec_prob;

        mio::ContactMatrixGroup<ScalarType> contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, num_agegroups);
        contact_matrix[0] =
            mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(num_agegroups, num_agegroups, 4.));
        model_group1.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        mio::SmootherCosine<ScalarType> smoothcos_prob(1.0);
        mio::StateAgeFunctionWrapper<ScalarType> prob(smoothcos_prob);
        model_group1.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[mio::AgeGroup(group)] = prob;
        model_group1.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[mio::AgeGroup(group)]   = prob;
        model_group1.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[mio::AgeGroup(group)]   = prob;
    }
    model_group1.parameters.set<mio::isecir::Seasonality>(0.);
    model_group1.parameters.set<mio::isecir::StartDay>(0);

    // Carry out simulation.
    mio::isecir::Simulation sim_group1(model_group1, dt);
    sim_group1.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated_group1    = sim_group1.get_result();
    mio::TimeSeries<ScalarType> transitions_simulated_group1 = sim_group1.get_transitions();

    /*****************************************************************************************************************/

    // Define vectors for compartments and transitions with expected results from example
    // (calculated by hand, see https://doi.org/10.1016/j.amc.2025.129636 for the used formulas).
    std::vector<ScalarType> secihurd_t0_baseline    = {4995., 0.5, 0., 4., 0., 0., 4990.5, 10.};
    std::vector<ScalarType> secihurd_t1_baseline    = {4994.00020016, 0.49989992, 0.49994996,    0.12498749,
                                                       1.03124687,    0.25781172, 4993.45699802, 10.12890586};
    std::vector<ScalarType> transitions_t1_baseline = {0.99979984, 0.99989992, 0.24997498, 0.24997498, 2.06249374,
                                                       2.06249374, 0.51562344, 0.51562344, 0.12890586, 0.12890586};

    // Define vectors with expected results for first model.
    int num_compartments = (int)mio::isecir::InfectionState::Count;
    std::vector<ScalarType> secihurd_t0_group0(num_compartments * num_agegroups, 0.);
    std::vector<ScalarType> secihurd_t1_group0(num_compartments * num_agegroups, 0.);
    std::vector<ScalarType> transitions_t1_group0(num_transitions * num_agegroups, 0.);

    // Copy expected results for first age group.
    std::copy(secihurd_t0_baseline.begin(), secihurd_t0_baseline.end(), secihurd_t0_group0.begin());
    std::copy(secihurd_t1_baseline.begin(), secihurd_t1_baseline.end(), secihurd_t1_group0.begin());
    std::copy(transitions_t1_baseline.begin(), transitions_t1_baseline.end(), transitions_t1_group0.begin());

    // Define vectors with expected results for second model.
    std::vector<ScalarType> secihurd_t0_group1(num_compartments * num_agegroups, 0.);
    std::vector<ScalarType> secihurd_t1_group1(num_compartments * num_agegroups, 0.);
    std::vector<ScalarType> transitions_t1_group1(num_transitions * num_agegroups, 0.);

    // Scale baseline vectors by baseline_scaling.
    std::transform(secihurd_t0_baseline.begin(), secihurd_t0_baseline.end(), secihurd_t0_baseline.begin(),
                   [baseline_scaling](ScalarType x) {
                       return baseline_scaling * x;
                   });
    std::transform(secihurd_t1_baseline.begin(), secihurd_t1_baseline.end(), secihurd_t1_baseline.begin(),
                   [baseline_scaling](ScalarType x) {
                       return baseline_scaling * x;
                   });
    std::transform(transitions_t1_baseline.begin(), transitions_t1_baseline.end(), transitions_t1_baseline.begin(),
                   [baseline_scaling](ScalarType x) {
                       return baseline_scaling * x;
                   });

    // Copy expected results for second age group.
    std::copy(secihurd_t0_baseline.begin(), secihurd_t0_baseline.end(), secihurd_t0_group1.begin() + num_compartments);
    std::copy(secihurd_t1_baseline.begin(), secihurd_t1_baseline.end(), secihurd_t1_group1.begin() + num_compartments);
    std::copy(transitions_t1_baseline.begin(), transitions_t1_baseline.end(),
              transitions_t1_group1.begin() + num_transitions);

    // Compare simulated compartments at time points t0 and t1.
    for (Eigen::Index i = 0; i < (Eigen::Index)secihurd_simulated_group0.get_num_elements(); i++) {
        EXPECT_NEAR(secihurd_simulated_group0[0][i], secihurd_t0_group0[i], 1e-8);
        EXPECT_NEAR(secihurd_simulated_group0[1][i], secihurd_t1_group0[i], 1e-8);
        EXPECT_NEAR(secihurd_simulated_group1[0][i], secihurd_t0_group1[i], 1e-8);
        EXPECT_NEAR(secihurd_simulated_group1[1][i], secihurd_t1_group1[i], 1e-8);
    }

    // Compare simulated transitions with expected results at time point t1.
    for (Eigen::Index i = 0; i < transitions_simulated_group0.get_num_elements(); i++) {
        EXPECT_NEAR(transitions_simulated_group0[transitions_simulated_group0.get_num_time_points() - 1][i],
                    transitions_t1_group0[i], 1e-8);
        EXPECT_NEAR(transitions_simulated_group1[transitions_simulated_group1.get_num_time_points() - 1][i],
                    transitions_t1_group1[i], 1e-8);
    }
}
