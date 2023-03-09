/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "load_test_data.h"
#include "ide_secir/model.h"
#include "ide_secir/parameters.h"
#include "ide_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>
#include <gtest/gtest.h>

class ModelTest : public testing::Test
{
protected:
    virtual void SetUp()
    {
        using Vec = mio::TimeSeries<ScalarType>::Vector;

        //Set initial conditions
        ScalarType N           = 10000;
        ScalarType Dead_before = 12;
        ScalarType dt          = 1;

        int num_transitions = (int)mio::isecir::InfectionTransitions::Count;

        Vec vec_init(num_transitions);
        mio::TimeSeries<ScalarType> init(num_transitions);
        vec_init << 25.0, 15.0, 8.0, 4.0, 1.0, 4.0, 1.0, 1.0, 1.0, 1.0;
        init.add_time_point(-10.0, vec_init);
        while (init.get_last_time() < 0) {
            vec_init *= 1.01;
            init.add_time_point(init.get_last_time() + dt, vec_init);
        }

        // Initialize model
        model = new mio::isecir::Model(std::move(init), 1, N, Dead_before);

        // Set working parameters.
        model->parameters.set<mio::isecir::TransitionDistributions>(
            std::vector<mio::isecir::DelayDistribution>(num_transitions, mio::isecir::DelayDistribution()));
        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransitions::Count, 0.5);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransitions::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransitions::ExposedToInfectedNoSymptoms)] = 1;
        model->parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);
        mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
        contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
        model->parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);
        model->parameters.set<mio::isecir::TransmissionProbabilityOnContact>(0.5);
        model->parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(0.5);
        model->parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(0.5);
    }

    virtual void TearDown()
    {
        delete model;
    }

public:
    mio::isecir::Model* model = nullptr;
};

// check if population stays constant over course of simulation
TEST_F(ModelTest, checkPopulationConservation)
{
    mio::TimeSeries<ScalarType> compartments = model->simulate(15);

    ScalarType num_persons_before = 0.0;
    ScalarType num_persons_after  = 0.0;

    for (auto i = 0; i < compartments[0].size(); i++) {
        num_persons_before += compartments[0][i];
        num_persons_after += compartments.get_last_value()[i];
    }

    EXPECT_NEAR(num_persons_after, num_persons_before, 1e-10);
}

// compare compartments with previous run
TEST_F(ModelTest, compareWithPreviousRun)
{
    auto compare                             = load_test_data_csv<ScalarType>("ide-secir-compare.csv");
    mio::TimeSeries<ScalarType> compartments = model->simulate(5);

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
TEST_F(ModelTest, compareWithPreviousRunTransitions)
{
    auto compare = load_test_data_csv<ScalarType>("ide-secir-transitions-compare.csv");
    model->simulate(5);
    auto transitions = model->get_transitions();

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

// check rsults of our simulation with example calculated by hand
// for exanḿple see Overleaf document
TEST(IdeSecir, checksimulationFunctions)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax        = 1;
    ScalarType N           = 10000;
    ScalarType Dead_before = 10;
    ScalarType dt          = 1;

    int num_transitions = (int)mio::isecir::InfectionTransitions::Count;

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    mio::TimeSeries<ScalarType> init(num_transitions);

    // add time points for initialization for transitions and death
    Vec vec_init(num_transitions);
    vec_init << 1.0, 0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    // add initial time point to time series
    init.add_time_point(-1, vec_init);
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), dt, N, Dead_before);

    // Set working parameters.
    model.parameters.set<mio::isecir::TransitionDistributions>(
        std::vector<mio::isecir::DelayDistribution>(num_transitions, mio::isecir::DelayDistribution()));
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransitions::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransitions::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransitions::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);
    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 1.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(0.5);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(1.0);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(1.0);

    // Carry out simulation.
    mio::TimeSeries<ScalarType> secihurd_simulated    = model.simulate(tmax);
    mio::TimeSeries<ScalarType> transitions_simulated = model.get_transitions();

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