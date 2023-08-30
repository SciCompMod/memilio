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
#include "memilio/epidemiology/state_age_function.h"
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
        mio::SmootherCosine smoothcos(2.0);
        mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
        std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);

        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        model->parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

        mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
        contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
        model->parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        mio::ExponentialDecay expdecay(0.5);
        mio::StateAgeFunctionWrapper prob(expdecay);
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
TEST(IdeSecir, checkSimulationFunctions)
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
    // In our example we use m_support_max = 2 for all DelayDistribution%s
    // For all TransitionDistribution%s we use a SmootherCosine Function with funcparam=2.
    // In this case, funcparam is equal to the support_max.
    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 2.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::SmootherCosine smoothcos_prob(2.0);
    mio::StateAgeFunctionWrapper prob(smoothcos_prob);
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

// a) Test if check_constraints() function correctly reports wrongly set parameters.
// b) Test if check_constraints() does not complain if parameters are set within correct ranges.
TEST(IdeSecir, testValueConstraints)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType N           = 10000;
    ScalarType Dead_before = 10;
    ScalarType dt          = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_transitions);

    // Add time points for initialization for transitions.
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
    // Add initial time point to time series.
    init.add_time_point(-12, vec_init);
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize a model.
    mio::isecir::Model model(std::move(init), N, Dead_before);

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Set wrong parameters and use check constraints.

    // Create invalid and valid function and first set invalid function for all parameters.
    // Go same order as in check_constraints().
    mio::ConstantFunction constant_func_neg(-1);
    mio::StateAgeFunctionWrapper prob_neg(constant_func_neg);
    mio::ConstantFunction constant_func_pos(1);
    mio::StateAgeFunctionWrapper prob_pos(constant_func_pos);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob_neg);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob_neg);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob_neg);

    // Warn, i.e., return true for wrong TransmissionProbabilityOnContact.
    auto constraint_check = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Correct wrong parameter so that next check can go through.
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob_pos);
    // Warn, i.e., return true for wrong RelativeTransmissionNoSymptoms.
    constraint_check = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Correct wrong parameter so that next check can go through.
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob_pos);
    // Warn, i.e., return true for wrong RiskOfInfectionFromSymptomatic.
    constraint_check = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Correct wrong parameter so that next check can go through.
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob_pos);

    // Set wrong values for InfectionTransitions, one after the other.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]          = 0.2;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)]   = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)]     = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)]   = 0.4;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)]        = -0.6;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);
    // Check all, stop at InfectedCriticalToDead.
    constraint_check = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check SusceptibleToExposed.
    model.parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        0.6;
    constraint_check = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check ExposedToInfectedNoSymptoms.
    model.parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::SusceptibleToExposed] = 1.0;
    model.parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = -1.0;
    constraint_check                                                       = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check sum InfectedNoSymptomsToInfectedSymptoms + InfectedNoSymptomsToRecovered.
    model.parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]   = 1.0;
    model.parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] = 0.9;
    constraint_check                                                         = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check sum InfectedSymptomsToInfectedSevere + InfectedSymptomsToRecovered.
    model.parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]    = 0.0;
    model.parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] = 0.2;
    constraint_check                                                            = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check sum InfectedSevereToInfectedCritical + InfectedCriticalToRecovered.
    model.parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] = 1.0;
    model.parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
        0.4;
    constraint_check = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check sum InfectedCriticalToDead + InfectedSevereToRecovered.
    model.parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
        0.0;
    model.parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        0.59;
    constraint_check = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Set wrong function type with unlimited support.
    model.parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        0.6;
    mio::ConstantFunction const_func(1.0);
    mio::StateAgeFunctionWrapper delaydistribution(const_func);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);
    constraint_check = model.parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check all parameter are correct.
    mio::ExponentialDecay expdecay(4.0);
    mio::StateAgeFunctionWrapper delaydistribution2(expdecay);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib2(num_transitions, delaydistribution2);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib2);
    constraint_check = model.parameters.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

// The idea of this test is to check whether the proportion between Recovered and Dead is as expected
// (after simulation for a long enough time, i.e. when the equlibrium is reached).
TEST(IdeSecir, checkProportionRecoveredDeath)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax        = 30;
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

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, Dead_before);

    // Set working parameters.
    // All TransitionDistribution%s are ExponentialDecay functions.
    // For all transitions we have funcparam=2 except for InfectedCriticalToRecovered where we set funcparam=3.
    mio::ExponentialDecay expdecay(4.0);
    mio::StateAgeFunctionWrapper delaydistribution(expdecay);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_parameter(3.0);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Set probabilities so that all individuals go from Susceptible to InfectedCritical with probability 1, from there they move
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

    mio::ExponentialDecay expdecay2(0.5);
    mio::StateAgeFunctionWrapper prob(expdecay2);
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
// (after simulation for a long enough time) does not change if we have a different m_support_max
// for the TransitionDistribution describing the transition from InfectedCritical To Recovered.
// We also check whether the equilibirum is reached earlier if m_support_max is chosen smaller.
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
    // Here we set the support_max for the TransitionDistribution%s differently for both models.
    // In both models, all TransitionDistribution%s are SmootherCosineFunctions

    // For model
    // All TransitionDistribution%s have funcparam=2
    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // For model2
    // All TransitionDistribution%s have funcparam=2 except fpr InfectedCriticalToRecovered wehre we set funcparam=7
    mio::SmootherCosine smoothcos2(2.0);
    mio::StateAgeFunctionWrapper delaydistribution2(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib2(num_transitions, delaydistribution2);
    vec_delaydistrib2[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_parameter(7.0);
    model2.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib2);

    // All remaining parameters are equal for both models.
    // Set probabilities so that all individuals go from Susceptible to InfectedCritical with probability 1, from there they move
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

    mio::ExponentialDecay expdecay(0.5);
    mio::StateAgeFunctionWrapper prob(expdecay);
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

    // Check whether equilibrium has been reached, only then it makes sense to compare results and times when equilibrium was reached.
    EXPECT_TRUE(secihurd_simulated[Eigen::Index(tmax / dt - 1)] == secihurd_simulated[Eigen::Index(tmax / dt - 2)]);
    EXPECT_TRUE(secihurd_simulated2[Eigen::Index(tmax / dt - 1)] == secihurd_simulated2[Eigen::Index(tmax / dt - 2)]);

    // Check whether both models have the same result at time tmax
    for (Eigen::Index i = 0; i < (Eigen::Index)mio::isecir::InfectionState::Count; i++) {
        EXPECT_NEAR(secihurd_simulated.get_last_value()[i], secihurd_simulated2.get_last_value()[i], 1e-8);
    }

    // Compute at what time the equilibrium was reached and check whether that time point is smaller for model than for model2
    // (as we have a smaller support_max in model than in model2)
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

TEST(IdeSecir, checkInfectionTransitions)
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
