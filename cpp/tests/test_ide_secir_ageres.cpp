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
    compare_compartments << 484.3056557672, 15.7685031055, 22.7020934123, 7.0615933479, 3.3491460693, 1.5803397070,
        4454.5548070034, 10.6778615873, 484.3056557672, 31.0934010790, 21.1271954388, 23.6370809253, 3.9106794140,
        1.7242153411, 4424.0110181177, 10.1907539167, 605.3820697090, 60.1973290710, 23.8046231705, 16.6085494134,
        3.6307172673, 1.6536810707, 4278.2949856871, 10.4280446109;

    ASSERT_EQ(compare_compartments.size(), static_cast<size_t>(compartments.get_last_value().size()));

    for (int j = 0; j < compare_compartments.size(); j++) {
        ASSERT_NEAR(compartments.get_last_value()[j], compare_compartments[j], 1e-7);
    }

    // Compare transitions at last time point with results from a previous run that are given here.

    mio::TimeSeries<ScalarType> transitions = sim.get_transitions();
    Eigen::VectorX<ScalarType> compare_transitions(num_transitions * num_agegroups);
    compare_transitions << 31.5370062111, 30.6497959470, 14.1231866958, 14.7543908776, 6.6982921386, 6.6982921386,
        3.1606794140, 3.1606794140, 1.4742153411, 1.4742153411, 31.5370062111, 29.5087817552, 14.7543908776,
        14.1231866958, 6.3213588280, 6.3213588280, 2.9484306823, 2.9484306823, 1.3533839877, 1.3533839877,
        39.4212577639, 30.1092463410, 14.4459531059, 14.4459531059, 6.5114345346, 6.5114345346, 3.0573621415,
        3.0573621415, 1.4155355888, 1.4155355888;

    ASSERT_EQ(compare_transitions.size(), static_cast<size_t>(transitions.get_last_value().size()));

    for (int j = 0; j < compare_transitions.size(); j++) {
        ASSERT_NEAR(transitions.get_last_value()[j], compare_transitions[j], 1e-7);
    }
}
