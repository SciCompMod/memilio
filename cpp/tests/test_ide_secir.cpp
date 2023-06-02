/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke, Anna Wendler
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

#include "boost/fusion/functional/invocation/invoke.hpp"
#include "load_test_data.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>
#include <gtest/gtest.h>

class ModelTestIdeSecir : public testing::Test
{
protected:
    virtual void SetUp()
    {
        using Vec = mio::TimeSeries<ScalarType>::Vector;

        //Set initial conditions
        ScalarType N           = 10000;
        ScalarType Dead_before = 12;

        int num_transitions = (int)mio::isecir::InfectionTransition::Count;

        Vec vec_init(num_transitions);
        mio::TimeSeries<ScalarType> init(num_transitions);
        vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 25.0;
        vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 1.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 1.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 1.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 1.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 1.0;
        init.add_time_point(-10.0, vec_init);
        while (init.get_last_time() < 0) {
            vec_init *= 1.01;
            init.add_time_point(init.get_last_time() + dt, vec_init);
        }

        // Initialize model
        model = new mio::isecir::Model(std::move(init), N, Dead_before);

        // Set working parameters.
        model->parameters.set<mio::isecir::TransitionDistributions>(
            std::vector<mio::isecir::DelayDistribution>(num_transitions, mio::isecir::DelayDistribution()));

        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        model->parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

        mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
        contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
        model->parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        mio::isecir::StateAgeFunctionWrapper prob;
        mio::isecir::ExponentialDecay expdecay(0.5);
        prob.set_state_age_function(expdecay);
        model->parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
        model->parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
        model->parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);
    }

    virtual void TearDown()
    {
        delete model;
    }

public:
    mio::isecir::Model* model = nullptr;
    ScalarType dt             = 1;
};

// check if population stays constant over course of simulation
TEST_F(ModelTestIdeSecir, checkPopulationConservation)
{
    mio::TimeSeries<ScalarType> compartments = simulate(0, 15, dt, *model);

    ScalarType num_persons_before = 0.0;
    ScalarType num_persons_after  = 0.0;

    for (auto i = 0; i < compartments[0].size(); i++) {
        num_persons_before += compartments[0][i];
        num_persons_after += compartments.get_last_value()[i];
    }

    EXPECT_NEAR(num_persons_after, num_persons_before, 1e-10);
}

// compare compartments with previous run
TEST_F(ModelTestIdeSecir, compareWithPreviousRun)
{
    auto compare                             = load_test_data_csv<ScalarType>("ide-secir-compare.csv");
    mio::TimeSeries<ScalarType> compartments = simulate(0, 5, dt, *model);

    ASSERT_EQ(compare.size(), static_cast<size_t>(compartments.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(compartments.get_num_elements()) + 1) << "at row " << i;
        ASSERT_NEAR(compartments.get_time(i), compare[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            ASSERT_NEAR(compartments.get_value(i)[j - 1], compare[i][j], 1e-7) << " at row " << i;
        }
    }
}

// compare transitions with previous run
TEST_F(ModelTestIdeSecir, compareWithPreviousRunTransitions)
{
    auto compare = load_test_data_csv<ScalarType>("ide-secir-transitions-compare.csv");

    mio::isecir::Simulation sim(*model, 0, dt);
    sim.advance(5);

    auto transitions = sim.get_transitions();

    size_t iter_0 = 0;
    while (transitions.get_time(iter_0) < compare[0][0]) {
        iter_0++;
    }

    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(transitions.get_num_elements()) + 1) << "at row " << i;
        ASSERT_NEAR(transitions.get_time(i + iter_0), compare[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            ASSERT_NEAR(transitions.get_value(i + iter_0)[j - 1], compare[i][j], 1e-7) << " at row " << i;
        }
    }
}

// check results of our simulation with an example calculated by hand
// for example see Overleaf document
TEST(IdeSecir, checksimulationFunctions)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax        = 1;
    ScalarType N           = 10000;
    ScalarType Dead_before = 10;
    ScalarType dt          = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    mio::TimeSeries<ScalarType> init(num_transitions);

    // add time points for initialization for transitions and death
    Vec vec_init(num_transitions);
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 0.0;
    // add initial time point to time series
    init.add_time_point(-1, vec_init);
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, Dead_before);

    // Set working parameters.
    // In our example we use m_max_support = 2 for all DelayDistribution%s
    std::vector<ScalarType> vec_max_support((int)mio::isecir::InfectionTransition::Count, 2);
    std::vector<mio::isecir::DelayDistribution> vec_delaydistrib(num_transitions, mio::isecir::DelayDistribution());
    for (int i = 0; i < (int)mio::isecir::InfectionTransition::Count; i++) {
        vec_delaydistrib[i].set_max_support(vec_max_support[i]);
    }
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::isecir::StateAgeFunctionWrapper prob;
    mio::isecir::SmootherCosine smoothcos(2);
    prob.set_state_age_function(smoothcos);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);

    // Carry out simulation.
    mio::isecir::Simulation sim(model, 0, dt);
    sim.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated    = sim.get_result();
    mio::TimeSeries<ScalarType> transitions_simulated = sim.get_transitions();

    // Define vectors with values from example (calculated by hand, see Overleaf document)
    Vec secihurd0((int)mio::isecir::InfectionState::Count);
    Vec secihurd1((int)mio::isecir::InfectionState::Count);
    Vec transitions1(num_transitions);
    secihurd0 << 4995, 0.5, 0, 4, 0, 0, 4990.5, 10;
    secihurd1 << 4994.00020016, 0.49989992, 0.49994996, 0.12498749, 1.03124687, 0.25781172, 4993.45699802, 10.12890586;
    transitions1 << 0.99979984, 0.99989991, 0.24997498, 0.24997498, 2.06249374, 2.06249374, 0.51562344, 0.51562344,
        0.12890586, 0.12890586;

    // Compare SECIHURD compartments at times 0 and 1
    for (Eigen::Index i = 0; i < (Eigen::Index)mio::isecir::InfectionState::Count; i++) {
        EXPECT_NEAR(secihurd_simulated[0][i], secihurd0[i], 1e-8);
        EXPECT_NEAR(secihurd_simulated[1][i], secihurd1[i], 1e-8);
    }

    // Compare transitions at time 1
    for (Eigen::Index i = 0; i < num_transitions; i++) {
        EXPECT_NEAR(transitions_simulated[transitions_simulated.get_num_time_points() - 1][i], transitions1[i], 1e-8);
    }
}

TEST(IdeSecir, checkStateAgeFunctionWrapper)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType N           = 10000;
    ScalarType Dead_before = 10;
    ScalarType dt          = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    mio::TimeSeries<ScalarType> init(num_transitions);

    // add time points for initialization for transitions
    Vec vec_init(num_transitions);
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 0.0;
    // add initial time point to time series
    init.add_time_point(-12, vec_init);
    // add further time points
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, Dead_before);

    // Check that StateAgeFunctions are only considered equal if they are of the same derived
    // class and have the same funcparam
    mio::isecir::ExponentialDecay expdecay(0.5);
    mio::isecir::ExponentialDecay expdecay2(0.5);
    mio::isecir::ExponentialDecay expdecay3(1.0);
    mio::isecir::SmootherCosine smoothcos(0.5);

    EXPECT_TRUE(expdecay == expdecay2);
    EXPECT_FALSE(expdecay == expdecay3);
    EXPECT_FALSE(expdecay == smoothcos);

    // Check that it also holds when a StateAgeFunctionWrapper is set with the respective functions
    mio::isecir::StateAgeFunctionWrapper prob(expdecay);
    mio::isecir::StateAgeFunctionWrapper prob2(expdecay2);
    mio::isecir::StateAgeFunctionWrapper prob3(expdecay3);
    mio::isecir::StateAgeFunctionWrapper prob4(smoothcos);

    EXPECT_TRUE(prob == prob2);
    EXPECT_FALSE(prob == prob3);
    EXPECT_FALSE(prob == prob4);

    // Check that it also holds when a parameter is set with the respective StateAgeFunctionWrapper
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob2);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob3);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob4);

    EXPECT_TRUE(prob == model.parameters.get<mio::isecir::TransmissionProbabilityOnContact>());
    EXPECT_FALSE(prob == model.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>());
    EXPECT_FALSE(prob == model.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>());
}

// The idea of this test is to check whether the proportion between Recovered and Dead is as expected
// (after simulation for a long enough time, i.e. when the equlibrium is reached).
TEST(IdeSecir, checkProportionRecoveredDeath)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax        = 20;
    ScalarType N           = 10000;
    ScalarType Dead_before = 10;
    ScalarType dt          = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    mio::TimeSeries<ScalarType> init(num_transitions);

    // add time points for initialization for transitions
    Vec vec_init(num_transitions);
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 0.0;
    // add initial time point to time series
    init.add_time_point(-12, vec_init);
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize two models.
    mio::isecir::Model model(std::move(init), N, Dead_before);

    // Set working parameters.
    std::vector<ScalarType> vec_max_support((int)mio::isecir::InfectionTransition::Count, 2);
    vec_max_support[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] = 3;
    std::vector<mio::isecir::DelayDistribution> vec_delaydistrib(num_transitions, mio::isecir::DelayDistribution());
    for (int i = 0; i < (int)mio::isecir::InfectionTransition::Count; i++) {
        vec_delaydistrib[i].set_max_support(vec_max_support[i]);
    }
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Set probabilties so that all indiviudal go from Susceptible to InfectedCritical with probability 1, from there they move
    // to Recovered or Dead with probability 0.5, respectively.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)]   = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)]     = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)]   = 0.4;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)]        = 0.6;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 1.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::isecir::StateAgeFunctionWrapper prob;
    mio::isecir::ExponentialDecay expdecay(0.5);
    prob.set_state_age_function(expdecay);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);

    // Carry out simulation.
    mio::isecir::Simulation sim(model, 0, dt);
    sim.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated = sim.get_result();

    // Check whether equilibrium has been reached, only then can we expect the right proportion
    // between Recovered and Dead
    EXPECT_TRUE(secihurd_simulated[Eigen::Index(tmax / dt - 1)] == secihurd_simulated[Eigen::Index(tmax / dt - 2)]);

    // Check whether equilibrium has the right proportion between Recovered and Dead
    EXPECT_NEAR((vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] /
                 vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)]) *
                    (secihurd_simulated.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered] -
                     secihurd_simulated[0][(Eigen::Index)mio::isecir::InfectionState::Recovered]),
                secihurd_simulated.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Dead] -
                    secihurd_simulated[0][(Eigen::Index)mio::isecir::InfectionState::Dead],
                1e-8);
}

// The idea of this test is to confirm that the equilibrium of the compartments
// (after simulation for a long enough time) does not change if we have a different m_max_support
// for the DelayDistribution describing the transition from InfectedCritical To Recovered.
// We also check whether the equilibirum is reached earlier if m_max_support is chosen smaller.
TEST(IdeSecir, compareEquilibria)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax        = 20;
    ScalarType N           = 10000;
    ScalarType Dead_before = 10;
    ScalarType dt          = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    mio::TimeSeries<ScalarType> init(num_transitions);

    // add time points for initialization for transitions
    Vec vec_init(num_transitions);
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 0.0;
    // add initial time point to time series
    init.add_time_point(-12, vec_init);
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    mio::TimeSeries<ScalarType> init2(init);

    // Initialize two models.
    mio::isecir::Model model(std::move(init), N, Dead_before);
    mio::isecir::Model model2(std::move(init2), N, Dead_before);

    // Set working parameters.
    // Here we set the max_support for the DelayDistribution differently for both models
    // For model
    std::vector<ScalarType> vec_max_support((int)mio::isecir::InfectionTransition::Count, 2);
    vec_max_support[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] = 2;
    std::vector<mio::isecir::DelayDistribution> vec_delaydistrib(num_transitions, mio::isecir::DelayDistribution());
    for (int i = 0; i < (int)mio::isecir::InfectionTransition::Count; i++) {
        vec_delaydistrib[i].set_max_support(vec_max_support[i]);
    }
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // For model2
    std::vector<ScalarType> vec_max_support2((int)mio::isecir::InfectionTransition::Count, 2);
    vec_max_support2[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] = 7;
    std::vector<mio::isecir::DelayDistribution> vec_delaydistrib2(num_transitions, mio::isecir::DelayDistribution());
    for (int i = 0; i < (int)mio::isecir::InfectionTransition::Count; i++) {
        vec_delaydistrib2[i].set_max_support(vec_max_support2[i]);
    }
    model2.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib2);

    // All remaining parameters are equal for both models.
    // Set probabilties so that all indiviudal go from Susceptible to InfectedCritical with probability 1, from there they move
    // to Recovered or Dead with probability 0.5, respectively.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)]   = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)]     = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)]   = 0.5;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)]        = 0.5;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);
    model2.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix                = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                     = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 1.));
    model.parameters.get<mio::isecir::ContactPatterns>()  = mio::UncertainContactMatrix(contact_matrix);
    model2.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::isecir::StateAgeFunctionWrapper prob;
    mio::isecir::ExponentialDecay expdecay(0.5);
    prob.set_state_age_function(expdecay);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);

    model2.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model2.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model2.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);

    // Carry out simulation.
    mio::isecir::Simulation sim(model, 0, dt);
    sim.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated = sim.get_result();

    mio::isecir::Simulation sim2(model2, 0, dt);
    sim2.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated2 = sim2.get_result();

    // Check whether both models have the same result at time tmax
    for (Eigen::Index i = 0; i < (Eigen::Index)mio::isecir::InfectionState::Count; i++) {
        EXPECT_NEAR(secihurd_simulated.get_last_value()[i], secihurd_simulated2.get_last_value()[i], 1e-8);
    }

    // Compute at what time the equilibrium was reached and check whether that time point is smaller for model than for model2
    // (as we have a smaller max_support in model than in model2)
    ScalarType equilibrium_time{};
    ScalarType equilibrium_time2{};
    for (int t = 0; t < secihurd_simulated.get_num_time_points() - 1; t++) {
        if (secihurd_simulated[t] == secihurd_simulated[t + 1]) {
            equilibrium_time = t;
            break;
        }
    }
    for (int t = 0; t < secihurd_simulated.get_num_time_points() - 1; t++) {
        if (secihurd_simulated2[t] == secihurd_simulated2[t + 1]) {
            equilibrium_time2 = t;
            break;
        }
    }

    EXPECT_TRUE(equilibrium_time <= equilibrium_time2);
}

TEST(IdeSecir, infection_transitions)
{
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.size(), mio::isecir::InfectionTransitionsCount);

    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(0),
              std::make_pair(mio::isecir::InfectionState::Susceptible, mio::isecir::InfectionState::Exposed));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(1),
              std::make_pair(mio::isecir::InfectionState::Exposed, mio::isecir::InfectionState::InfectedNoSymptoms));
    EXPECT_EQ(
        mio::isecir::InfectionTransitionsMap.at(2),
        std::make_pair(mio::isecir::InfectionState::InfectedNoSymptoms, mio::isecir::InfectionState::InfectedSymptoms));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(3),
              std::make_pair(mio::isecir::InfectionState::InfectedNoSymptoms, mio::isecir::InfectionState::Recovered));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(4), std::make_pair(mio::isecir::InfectionState::InfectedSymptoms,
                                                                         mio::isecir::InfectionState::InfectedSevere));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(5),
              std::make_pair(mio::isecir::InfectionState::InfectedSymptoms, mio::isecir::InfectionState::Recovered));
    EXPECT_EQ(
        mio::isecir::InfectionTransitionsMap.at(6),
        std::make_pair(mio::isecir::InfectionState::InfectedSevere, mio::isecir::InfectionState::InfectedCritical));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(7),
              std::make_pair(mio::isecir::InfectionState::InfectedSevere, mio::isecir::InfectionState::Recovered));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(8),
              std::make_pair(mio::isecir::InfectionState::InfectedCritical, mio::isecir::InfectionState::Dead));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(9),
              std::make_pair(mio::isecir::InfectionState::InfectedCritical, mio::isecir::InfectionState::Recovered));
}