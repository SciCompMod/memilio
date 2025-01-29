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

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <gtest/gtest.h>
#include "smm/model.h"
#include "smm/parameters.h"
#include "smm/simulation.h"
#include "abm_helpers.h"
#include <iostream>
#include <numeric>
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

    //Set adoption rates
    std::vector<mio::smm::AdoptionRate<InfectionState>> adoption_rates;
    adoption_rates.push_back({InfectionState::S,
                              InfectionState::E,
                              mio::smm::Region(0),
                              0.1,
                              {InfectionState::C, InfectionState::I},
                              {1, 0.5}});
    adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::smm::Region(0), 0.2, {}, {}});

    //Initialize model populations
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
    //Same test as 'evaluateAdoptionRate' only for spatial transition rates
    using Model = mio::smm::Model<2, InfectionState>;

    Model model;
    //Initialize model populations
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
    //Set transition rates
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
    //Initialize model populations
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
    //Set adoption and transition rates
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

    //Mock exponential distribution to control the normalized waiting times that are drawn
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::AtLeast(7))
        .WillOnce(Return(0.5)) //adoption event S->E
        .WillOnce(Return(0.5)) //E->C
        .WillOnce(Return(0.5)) //C->R
        .WillOnce(Return(0.5)) //C->I
        .WillOnce(Return(0.0397)) //I->R
        .WillOnce(Return(0.5)) //I->D
        .WillOnce(Return(0.0031)) //spatial transition event 1->0
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
    //adoption from I to R in second time step
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::I)], 0);
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::R)], 1);
    EXPECT_EQ(sim.get_result().get_value(
                  1)[static_cast<size_t>(InfectionState::Count) + static_cast<size_t>(InfectionState::R)],
              1);
    //spatial transition in third time step
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::I)], 0);
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::R)], 2);
}

TEST(TestSMMSimulation, stopsAtTmax)
{
    using testing::Return;
    using Model = mio::smm::Model<2, InfectionState>;

    Model model;

    //Set adoption and spatial transition rates
    std::vector<mio::smm::AdoptionRate<InfectionState>> adoption_rates;
    std::vector<mio::smm::TransitionRate<InfectionState>> transition_rates;

    adoption_rates.push_back({InfectionState::S,
                              InfectionState::E,
                              mio::smm::Region(0),
                              0.1,
                              {InfectionState::C, InfectionState::I},
                              {1, 0.5}});

    transition_rates.push_back({InfectionState::R, mio::smm::Region(1), mio::smm::Region(0), 0.01});

    // model.parameters.get<mio::smm::AdoptionRates<InfectionState>>()   = adoption_rates;
    // model.parameters.get<mio::smm::TransitionRates<InfectionState>>() = transition_rates;

    //Simulate for 30 days
    // auto sim = mio::Simulation<double, Model>(model, 0.0, 0.1);
    // sim.advance(30.);
    // //As model populations are all zero only t0 and tmax should be logged
    // EXPECT_EQ(sim.get_result().get_num_time_points(), 2);
    // EXPECT_EQ(sim.get_result().get_last_time(), 30.);
}

TEST(TestSMMSimulation, covergence)
{
    using testing::Return;
    using Model = mio::smm::Model<2, InfectionState>;

    double pop1      = 1000;
    double pop2      = 10000;
    size_t num_runs1 = 100;
    size_t num_runs2 = 10000;
    double rate      = 0.01;

    //Only set one spatial transition rate
    std::vector<mio::smm::TransitionRate<InfectionState>> transition_rates(
        1, {InfectionState::S, mio::smm::Region(0), mio::smm::Region(1), rate});

    std::vector<double> transitions1(num_runs1);
    std::vector<double> transitions2(num_runs2);

    //First try 100 unit-time step simulation with 1000 agents
    for (size_t n = 0; n < num_runs1; ++n) {
        Model model;

        model.populations[{mio::smm::Region(0), InfectionState::S}] = pop1;
        model.populations[{mio::smm::Region(0), InfectionState::E}] = 0;
        model.populations[{mio::smm::Region(0), InfectionState::C}] = 0;
        model.populations[{mio::smm::Region(0), InfectionState::I}] = 0;
        model.populations[{mio::smm::Region(0), InfectionState::R}] = 0;
        model.populations[{mio::smm::Region(0), InfectionState::D}] = 0;

        model.populations[{mio::smm::Region(1), InfectionState::S}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::E}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::C}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::I}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::R}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::D}]       = 0;
        model.parameters.get<mio::smm::TransitionRates<InfectionState>>() = transition_rates;
        model.get_rng().seed({static_cast<uint32_t>(n)});
        auto sim = mio::Simulation<double, Model>(model, 0.0, 1.0);
        sim.advance(1.);
        transitions1[n] = sim.get_model().populations[{mio::smm::Region(1), InfectionState::S}];
    }

    //Then try 10000 unit-time step simulation with 10000 agents
    for (size_t n = 0; n < num_runs2; ++n) {
        Model model;

        model.populations[{mio::smm::Region(0), InfectionState::S}] = pop2;
        model.populations[{mio::smm::Region(0), InfectionState::E}] = 0;
        model.populations[{mio::smm::Region(0), InfectionState::C}] = 0;
        model.populations[{mio::smm::Region(0), InfectionState::I}] = 0;
        model.populations[{mio::smm::Region(0), InfectionState::R}] = 0;
        model.populations[{mio::smm::Region(0), InfectionState::D}] = 0;

        model.populations[{mio::smm::Region(1), InfectionState::S}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::E}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::C}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::I}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::R}]       = 0;
        model.populations[{mio::smm::Region(1), InfectionState::D}]       = 0;
        model.parameters.get<mio::smm::TransitionRates<InfectionState>>() = transition_rates;
        model.get_rng().seed({static_cast<uint32_t>(n)});
        auto sim = mio::Simulation<double, Model>(model, 0.0, 1.0);
        sim.advance(1.);
        transitions2[n] = sim.get_model().populations[{mio::smm::Region(1), InfectionState::S}];
    }
    //The number of transitions from region 0 to region 1 should be approx. rate * pop in region 0
    double rel_diff1 =
        std::abs(rate * pop1 - std::abs(std::accumulate(transitions1.begin(), transitions1.end(), 0.0) / num_runs1)) /
        (rate * pop1);
    double rel_diff2 =
        std::abs(rate * pop2 - std::abs(std::accumulate(transitions2.begin(), transitions2.end(), 0.0) / num_runs2)) /
        (rate * pop2);

    EXPECT_GE(rel_diff1, rel_diff2);
}
