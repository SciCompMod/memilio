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
//#include "ide_secir/model.cpp"
#include "ide_secir/parameters.h"
#include "ide_secir/infection_state.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/uncertain_matrix.h"
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
        init.add_time_point(-10, vec_init);
        while (init.get_last_time() < 0) {
            init.add_time_point(init.get_last_time() + dt, init.get_last_value() * 1.01);
        }


        model = new mio::isecir::Model(std::move(init), dt, N, Dead0);

        // Set working parameters.
        // model->parameters.set<mio::isecir::TransitionDistributions>(
        //     std::vector<mio::isecir::DelayDistribution>(num_transitions, mio::isecir::DelayDistribution()));
        // model->parameters.set<mio::isecir::TransitionProbabilities>(std::vector<ScalarType>(num_transitions, 0.5));
        // mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
        // contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
        // model->parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);
        // model->parameters.set<mio::isecir::TransmissionProbabilityOnContact>(1.0);
        // model->parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(1.0);
        // model->parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(1.0);
    }

    virtual void TearDown()
    {
        delete model;
    }

public:
    mio::isecir::Model* model = nullptr;
};

TEST_F(ModelTest, checkPopulationConservation)
{
    mio::TimeSeries<ScalarType> compartments = model->simulate(15);

    ScalarType num_persons_0 = 0.0;

    for (auto i = 0; i < compartments[0].size(); i++) {
        num_persons_0 += compartments[0][i];
    }

    ScalarType num_persons = 0.0;

    for (auto i = 0; i < compartments.get_last_value().size(); i++) {
        num_persons += compartments.get_last_value()[i];
    }

    EXPECT_NEAR(num_persons,num_persons_0,1e-10);
    // using Vec = mio::TimeSeries<ScalarType>::Vector;

    // ScalarType tmax  = 10;
    // ScalarType N     = 10000;
    // ScalarType Dead0 = 12;
    // ScalarType dt    = 1;

    // int num_transitions = (int)mio::isecir::InfectionTransitions::Count;

    // create TimeSeries with num_transitions elements where transitions needed for simulation will be stored
    // mio::TimeSeries<ScalarType> init(num_transitions);

    // add time points for initialization
    // Vec vec_init(num_transitions);
    // vec_init << 30.0, 15.0, 8.0, 4.0, 1.0, 4.0, 1.0, 1.0, 1.0, 1.0;
    // init.add_time_point(-10, vec_init);
    // while (init.get_last_time() < 0) {
    //     init.add_time_point(init.get_last_time() + dt, init.get_last_value() * 1.01);
    // }

    // Initialize model.
    // mio::isecir::Model model(std::move(init), dt, N, Dead0);
    // ASSERT_EQ(1, 1);
    //model->simulate(15);
    // auto sim_result = model->calculate_EIR();
}

// TEST(IdeSecir, testParamConstructors)
// {
// }

// TEST(IdeSecir, infection_transitions)
// {
//     EXPECT_EQ(mio::isecir::InfectionTransitionsMap.size(), mio::isecir::InfectionTransitionsCount);

//     EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(0),
//               std::make_pair(mio::isecir::InfectionState::Susceptible, mio::isecir::InfectionState::Exposed));
//     EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(1),
//               std::make_pair(mio::isecir::InfectionState::Exposed, mio::isecir::InfectionState::InfectedNoSymptoms));
//     EXPECT_EQ(
//         mio::isecir::InfectionTransitionsMap.at(2),
//         std::make_pair(mio::isecir::InfectionState::InfectedNoSymptoms, mio::isecir::InfectionState::InfectedSymptoms));
//     EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(3),
//               std::make_pair(mio::isecir::InfectionState::InfectedNoSymptoms, mio::isecir::InfectionState::Recovered));
//     EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(4), std::make_pair(mio::isecir::InfectionState::InfectedSymptoms,
//                                                                          mio::isecir::InfectionState::InfectedSevere));
//     EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(5),
//               std::make_pair(mio::isecir::InfectionState::InfectedSymptoms, mio::isecir::InfectionState::Recovered));
//     EXPECT_EQ(
//         mio::isecir::InfectionTransitionsMap.at(6),
//         std::make_pair(mio::isecir::InfectionState::InfectedSevere, mio::isecir::InfectionState::InfectedCritical));
//     EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(7),
//               std::make_pair(mio::isecir::InfectionState::InfectedSevere, mio::isecir::InfectionState::Recovered));
//     EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(8),
//               std::make_pair(mio::isecir::InfectionState::InfectedCritical, mio::isecir::InfectionState::Dead));
//     EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(9),
//               std::make_pair(mio::isecir::InfectionState::InfectedCritical, mio::isecir::InfectionState::Recovered));
// }