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

//#include "load_test_data.h"
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

        ScalarType N     = 10000;
        ScalarType Dead0 = 12;
        ScalarType dt    = 1;

        int num_transitions = (int)mio::isecir::InfectionTransitions::Count;

        Vec vec_init(num_transitions);
        mio::TimeSeries<ScalarType> init(num_transitions);
        vec_init << 30.0, 15.0, 8.0, 4.0, 1.0, 4.0, 1.0, 1.0, 1.0, 1.0;
        init.add_time_point(-10.0, vec_init);
        while (init.get_last_time() < 0) {
            vec_init *=1.01;
            init.add_time_point(init.get_last_time() + dt, vec_init);
        }

        model = new mio::isecir::Model(std::move(init), 1, N, Dead0);

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

TEST_F(ModelTest, printTransitions)
{
    std::cout << " In TEST_F: \n";
    model->print_transitions();
    std::cout << "\n" << std::endl;
}

TEST_F(ModelTest, checkPopulationConservation)
{
    mio::TimeSeries<ScalarType> compartments = model->simulate(15);

    ScalarType num_persons_before = 0.0;
    ScalarType num_persons_after = 0.0;

    for (auto i = 0; i < compartments[0].size(); i++) {
        num_persons_before += compartments[0][i];
        num_persons_after += compartments.get_last_value()[i];
    }

    EXPECT_NEAR(num_persons_after,num_persons_before,1e-10);
}

TEST_F(ModelTest, checkPrivateSimulationFunctions)
{
    // compute_susceptibles

    // update_forceofinfection

    // compute_flow

    // flows_current_timestep

    // compute_totaldeaths

    // compute_recovered

    // compute_compartment

    // compartments_current_timestep_ECIHU
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