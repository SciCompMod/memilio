/* 
* Copyright (C) 2020-2024 German Aerospace Center (DLR-SC)
*
* Authors: Julia Bicker
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

#include <cstddef>
#include <gtest/gtest.h>
#include "smm/model.h"
#include "smm/parameters.h"
#include "smm/simulation.h"
#include "abm_helpers.h"
#include <vector>

enum class InfectionState
{
    S,
    E,
    C,
    I,
    R,
    D,
    Count

};

TEST(TestSMM, evaluateAdoptionRate)
{
    using Model = mio::smm::Model<1, InfectionState>;

    Model model;

    std::vector<mio::smm::AdoptionRate<InfectionState>> adoption_rates;
    adoption_rates.push_back({InfectionState::S,
                              InfectionState::E,
                              mio::smm::Region(0),
                              0.1,
                              {InfectionState::C, InfectionState::I},
                              {1, 0.5}});
    adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::smm::Region(0), 0.2, {}, {}});

    model.populations[{mio::smm::Region(0), InfectionState::S}] = 50;
    model.populations[{mio::smm::Region(0), InfectionState::E}] = 10;
    model.populations[{mio::smm::Region(0), InfectionState::C}] = 5;
    model.populations[{mio::smm::Region(0), InfectionState::I}] = 0;
    model.populations[{mio::smm::Region(0), InfectionState::R}] = 0;
    model.populations[{mio::smm::Region(0), InfectionState::D}] = 0;

    EXPECT_EQ(model.evaluate(adoption_rates[0], model.populations.get_compartments()), 5. / 13.);
    EXPECT_EQ(model.evaluate(adoption_rates[1], model.populations.get_compartments()), 2.);
}

TEST(TestSMM, evaluateTransitionRate)
{
    using Model = mio::smm::Model<2, InfectionState>;

    Model model;

    model.populations[{mio::smm::Region(0), InfectionState::S}] = 50;
    model.populations[{mio::smm::Region(0), InfectionState::E}] = 10;
    model.populations[{mio::smm::Region(0), InfectionState::C}] = 5;
    model.populations[{mio::smm::Region(0), InfectionState::I}] = 0;
    model.populations[{mio::smm::Region(0), InfectionState::R}] = 0;
    model.populations[{mio::smm::Region(0), InfectionState::D}] = 0;

    model.populations[{mio::smm::Region(1), InfectionState::S}] = 55;
    model.populations[{mio::smm::Region(1), InfectionState::E}] = 10;
    model.populations[{mio::smm::Region(1), InfectionState::C}] = 0;
    model.populations[{mio::smm::Region(1), InfectionState::I}] = 0;
    model.populations[{mio::smm::Region(1), InfectionState::R}] = 0;
    model.populations[{mio::smm::Region(1), InfectionState::D}] = 0;

    std::vector<mio::smm::TransitionRate<InfectionState>> transition_rates;
    transition_rates.push_back({InfectionState::S, mio::smm::Region(0), mio::smm::Region(1), 0.01});
    transition_rates.push_back({InfectionState::E, mio::smm::Region(1), mio::smm::Region(0), 0.1});

    EXPECT_EQ(model.evaluate(transition_rates[0], model.populations.get_compartments()), 0.5);
    EXPECT_EQ(model.evaluate(transition_rates[1], model.populations.get_compartments()), 1.);
}

TEST(TestSMMSimulation, advance)
{
    using testing::Return;
    using Model = mio::smm::Model<2, InfectionState>;

    Model model;

    model.populations[{mio::smm::Region(0), InfectionState::S}] = 1;
    model.populations[{mio::smm::Region(0), InfectionState::E}] = 0;
    model.populations[{mio::smm::Region(0), InfectionState::C}] = 0;
    model.populations[{mio::smm::Region(0), InfectionState::I}] = 1;
    model.populations[{mio::smm::Region(0), InfectionState::R}] = 0;
    model.populations[{mio::smm::Region(0), InfectionState::D}] = 0;

    model.populations[{mio::smm::Region(1), InfectionState::S}] = 0;
    model.populations[{mio::smm::Region(1), InfectionState::E}] = 0;
    model.populations[{mio::smm::Region(1), InfectionState::C}] = 0;
    model.populations[{mio::smm::Region(1), InfectionState::I}] = 0;
    model.populations[{mio::smm::Region(1), InfectionState::R}] = 1;
    model.populations[{mio::smm::Region(1), InfectionState::D}] = 0;

    std::vector<mio::smm::AdoptionRate<InfectionState>> adoption_rates;
    std::vector<mio::smm::TransitionRate<InfectionState>> transition_rates;

    adoption_rates.push_back({InfectionState::S,
                              InfectionState::E,
                              mio::smm::Region(0),
                              0.1,
                              {InfectionState::C, InfectionState::I},
                              {1, 0.5}});
    adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::smm::Region(0), 1.0 / 5., {}, {}});
    adoption_rates.push_back({InfectionState::C, InfectionState::R, mio::smm::Region(0), 0.2 / 3., {}, {}});
    adoption_rates.push_back({InfectionState::C, InfectionState::I, mio::smm::Region(0), 0.8 / 3., {}, {}});
    adoption_rates.push_back({InfectionState::I, InfectionState::R, mio::smm::Region(0), 0.99 / 5., {}, {}});
    adoption_rates.push_back({InfectionState::I, InfectionState::D, mio::smm::Region(0), 0.01 / 5., {}, {}});

    transition_rates.push_back({InfectionState::R, mio::smm::Region(1), mio::smm::Region(0), 0.01});

    model.parameters.get<mio::smm::AdoptionRates<InfectionState>>()   = adoption_rates;
    model.parameters.get<mio::smm::TransitionRates<InfectionState>>() = transition_rates;

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::AtLeast(7))
        .WillOnce(Return(0.5))
        .WillOnce(Return(0.5))
        .WillOnce(Return(0.5))
        .WillOnce(Return(0.5))
        .WillOnce(Return(0.0397))
        .WillOnce(Return(0.5))
        .WillOnce(Return(0.0031))
        .WillRepeatedly(testing::Return(1.0));

    auto sim = mio::Simulation<double, Model>(model, 0.0, 0.1);
    sim.advance(30.);
    //initial values
    EXPECT_EQ(sim.get_result().get_value(0)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(0)[static_cast<size_t>(InfectionState::I)], 1);
    EXPECT_EQ(sim.get_result().get_value(
                  0)[static_cast<size_t>(InfectionState::Count) + static_cast<size_t>(InfectionState::R)],
              1);
    //no event happens in first time step
    EXPECT_GE(sim.get_result().get_time(1), 0.2);
    //adoption from I to R
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::I)], 0);
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::R)], 1);
    EXPECT_EQ(sim.get_result().get_value(
                  1)[static_cast<size_t>(InfectionState::Count) + static_cast<size_t>(InfectionState::R)],
              1);
    //spatial transition
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::I)], 0);
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::R)], 2);
    EXPECT_EQ(sim.get_result().get_last_time(), 30.);
}
