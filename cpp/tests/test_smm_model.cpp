/*
* Copyright (C) 2020-2025 MEmilio
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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <gtest/gtest.h>
#include "memilio/utils/compiler_diagnostics.h"
#include "smm/model.h"
#include "smm/parameters.h"
#include "smm/simulation.h"
#include "abm_helpers.h"
#include "memilio/epidemiology/adoption_rate.h"
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
    //Test whether the adoption rates are evaluated correctly.
    //First-order adoption rates are given by:
    // rate.factor * N(rate.from, rate.region)
    // with N(from, region) the population in Region "region" having infection state "from"
    //Second-order adoption rates are given by:
    // rate.factor * N(rate.from, rate.region)/total_pop * sum (over all rate.influences) influence.factor * N(influence.status, rate.region)
    using Model = mio::smm::Model<double, 1, InfectionState>;

    Model model;

    //Set adoption rates
    std::vector<mio::AdoptionRate<double, InfectionState>> adoption_rates;
    //Second-order adoption
    adoption_rates.push_back({InfectionState::S,
                              InfectionState::E,
                              mio::regions::Region(0),
                              0.1,
                              {{InfectionState::C, 1}, {InfectionState::I, 0.5}}});
    //First-order adoption
    adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::regions::Region(0), 0.2, {}});

    //Initialize model populations
    model.populations[{mio::regions::Region(0), InfectionState::S}] = 50;
    model.populations[{mio::regions::Region(0), InfectionState::E}] = 10;
    model.populations[{mio::regions::Region(0), InfectionState::C}] = 5;
    model.populations[{mio::regions::Region(0), InfectionState::I}] = 0;
    model.populations[{mio::regions::Region(0), InfectionState::R}] = 0;
    model.populations[{mio::regions::Region(0), InfectionState::D}] = 0;

    EXPECT_EQ(model.evaluate(adoption_rates[0], model.populations.get_compartments()), 5. / 13.);
    EXPECT_EQ(model.evaluate(adoption_rates[1], model.populations.get_compartments()), 2.);
}

TEST(TestSMM, evaluateTransitionRate)
{
    //Same test as 'evaluateAdoptionRate' only for spatial transition rates.
    //Transition rates are given by: rate.factor * N(rate.status, rate.from)
    using Model = mio::smm::Model<double, 2, InfectionState>;

    Model model;
    //Initialize model populations
    model.populations[{mio::regions::Region(0), InfectionState::S}] = 50;
    model.populations[{mio::regions::Region(0), InfectionState::E}] = 10;
    model.populations[{mio::regions::Region(0), InfectionState::C}] = 5;
    model.populations[{mio::regions::Region(0), InfectionState::I}] = 0;
    model.populations[{mio::regions::Region(0), InfectionState::R}] = 0;
    model.populations[{mio::regions::Region(0), InfectionState::D}] = 0;

    model.populations[{mio::regions::Region(1), InfectionState::S}] = 55;
    model.populations[{mio::regions::Region(1), InfectionState::E}] = 10;
    model.populations[{mio::regions::Region(1), InfectionState::C}] = 0;
    model.populations[{mio::regions::Region(1), InfectionState::I}] = 0;
    model.populations[{mio::regions::Region(1), InfectionState::R}] = 0;
    model.populations[{mio::regions::Region(1), InfectionState::D}] = 0;
    //Set transition rates
    std::vector<mio::smm::TransitionRate<double, InfectionState>> transition_rates;
    transition_rates.push_back({InfectionState::S, mio::regions::Region(0), mio::regions::Region(1), 0.01});
    transition_rates.push_back({InfectionState::E, mio::regions::Region(1), mio::regions::Region(0), 0.1});

    EXPECT_EQ(model.evaluate(transition_rates[0], model.populations.get_compartments()), 0.5);
    EXPECT_EQ(model.evaluate(transition_rates[1], model.populations.get_compartments()), 1.);
}

TEST(TestSMMSimulation, advance)
{
    //Test whether Gillespie algorithm calculates events in the correct order
    using testing::Return;
    using Model = mio::smm::Model<double, 2, InfectionState>;

    Model model;
    //Initialize model populations
    model.populations[{mio::regions::Region(0), InfectionState::S}] = 1;
    model.populations[{mio::regions::Region(0), InfectionState::E}] = 0;
    model.populations[{mio::regions::Region(0), InfectionState::C}] = 0;
    model.populations[{mio::regions::Region(0), InfectionState::I}] = 1;
    model.populations[{mio::regions::Region(0), InfectionState::R}] = 0;
    model.populations[{mio::regions::Region(0), InfectionState::D}] = 0;

    model.populations[{mio::regions::Region(1), InfectionState::S}] = 0;
    model.populations[{mio::regions::Region(1), InfectionState::E}] = 0;
    model.populations[{mio::regions::Region(1), InfectionState::C}] = 0;
    model.populations[{mio::regions::Region(1), InfectionState::I}] = 0;
    model.populations[{mio::regions::Region(1), InfectionState::R}] = 1;
    model.populations[{mio::regions::Region(1), InfectionState::D}] = 0;
    //Set adoption and transition rates
    std::vector<mio::AdoptionRate<double, InfectionState>> adoption_rates;
    std::vector<mio::smm::TransitionRate<double, InfectionState>> transition_rates;

    //Second-order adoption
    adoption_rates.push_back({InfectionState::S,
                              InfectionState::E,
                              mio::regions::Region(0),
                              0.1,
                              {{InfectionState::C, 1}, {InfectionState::I, 0.5}}});
    //First-order adoptions
    adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::regions::Region(0), 1.0 / 5., {}});
    adoption_rates.push_back({InfectionState::C, InfectionState::R, mio::regions::Region(0), 0.2 / 3., {}});
    adoption_rates.push_back({InfectionState::C, InfectionState::I, mio::regions::Region(0), 0.8 / 3., {}});
    adoption_rates.push_back({InfectionState::I, InfectionState::R, mio::regions::Region(0), 0.99 / 5., {}});
    adoption_rates.push_back({InfectionState::I, InfectionState::D, mio::regions::Region(0), 0.01 / 5., {}});

    //Spatial transition
    transition_rates.push_back({InfectionState::R, mio::regions::Region(1), mio::regions::Region(0), 0.01});

    model.parameters.get<mio::smm::AdoptionRates<double, InfectionState>>()   = adoption_rates;
    model.parameters.get<mio::smm::TransitionRates<double, InfectionState>>() = transition_rates;

    //Mock exponential distribution to control the normalized waiting times that are drawn
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::AtLeast(7))
        .WillOnce(Return(0.5)) //adoption event S->E
        .WillOnce(Return(0.5)) //E->C
        .WillOnce(Return(0.5)) //C->R
        .WillOnce(Return(0.5)) //C->I
        .WillOnce(Return(0.0397)) //I->R, corresponds to global time 0.20050
        .WillOnce(Return(0.5)) //I->D
        .WillOnce(Return(0.0031)) //spatial transition event 1->0, corresponds to global time 0.31
        .WillRepeatedly(testing::Return(1.0));

    auto sim = mio::smm::Simulation(model, 0.0, 0.1);
    sim.advance(30.);
    //Check whether first value in result time series corresponds to initial values
    EXPECT_EQ(sim.get_result().get_value(0)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(0)[static_cast<size_t>(InfectionState::I)], 1);
    EXPECT_EQ(sim.get_result().get_value(
                  0)[static_cast<size_t>(InfectionState::Count) + static_cast<size_t>(InfectionState::R)],
              1);
    //no event happens in first two time steps i.e. the first saved result is after t=0.2
    EXPECT_GE(sim.get_result().get_time(1), 0.2);
    //adoption from I to R is first event
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::I)], 0);
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::R)], 1);
    EXPECT_EQ(sim.get_result().get_value(
                  1)[static_cast<size_t>(InfectionState::Count) + static_cast<size_t>(InfectionState::R)],
              1);
    //spatial transition is second event
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::I)], 0);
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::R)], 2);
}

TEST(TestSMMSimulation, stopsAtTmax)
{
    //Test whether simulation stops at tmax and whether system state at tmax is saved if there are no adoptions/transitions
    using testing::Return;
    using Model = mio::smm::Model<double, 2, InfectionState>;

    Model model;

    //Set adoption and spatial transition rates
    std::vector<mio::AdoptionRate<double, InfectionState>> adoption_rates;
    std::vector<mio::smm::TransitionRate<double, InfectionState>> transition_rates;

    adoption_rates.push_back({InfectionState::S,
                              InfectionState::E,
                              mio::regions::Region(0),
                              0.1,
                              {{InfectionState::C, 1}, {InfectionState::I, 0.5}}});

    transition_rates.push_back({InfectionState::R, mio::regions::Region(1), mio::regions::Region(0), 0.01});

    model.parameters.get<mio::smm::AdoptionRates<double, InfectionState>>()   = adoption_rates;
    model.parameters.get<mio::smm::TransitionRates<double, InfectionState>>() = transition_rates;

    //As populations are not set they have value 0 i.e. no events will happen
    //Simulate for 30 days
    auto sim = mio::smm::Simulation(model, 0.0, 0.1);
    sim.advance(30.);
    //As model populations are all zero only t0 and tmax should be logged
    EXPECT_EQ(sim.get_result().get_num_time_points(), 2);
    EXPECT_EQ(sim.get_result().get_last_time(), 30.);
}

TEST(TestSMMSimulation, convergence)
{
    //Test whether the mean number of transitions in one unit-time step corresponds to the expected number of transitions
    // given by rate * pop up to some tolerance
    using testing::Return;
    using Model = mio::smm::Model<double, 2, InfectionState>;

    double pop      = 1000;
    size_t num_runs = 100;
    double rate     = 0.2;

    //Only set one spatial transition rate
    std::vector<mio::smm::TransitionRate<double, InfectionState>> transition_rates(
        1, {InfectionState::S, mio::regions::Region(0), mio::regions::Region(1), rate});

    double expected_num_trans = rate * pop;
    std::vector<double> transitions(num_runs);

    //Do 100 unit-time step simulations with 1000 agents
    for (size_t n = 0; n < num_runs; ++n) {
        Model model;

        model.populations[{mio::regions::Region(0), InfectionState::S}] = pop;
        model.populations[{mio::regions::Region(0), InfectionState::E}] = 0;
        model.populations[{mio::regions::Region(0), InfectionState::C}] = 0;
        model.populations[{mio::regions::Region(0), InfectionState::I}] = 0;
        model.populations[{mio::regions::Region(0), InfectionState::R}] = 0;
        model.populations[{mio::regions::Region(0), InfectionState::D}] = 0;

        model.populations[{mio::regions::Region(1), InfectionState::S}]           = 0;
        model.populations[{mio::regions::Region(1), InfectionState::E}]           = 0;
        model.populations[{mio::regions::Region(1), InfectionState::C}]           = 0;
        model.populations[{mio::regions::Region(1), InfectionState::I}]           = 0;
        model.populations[{mio::regions::Region(1), InfectionState::R}]           = 0;
        model.populations[{mio::regions::Region(1), InfectionState::D}]           = 0;
        model.parameters.get<mio::smm::TransitionRates<double, InfectionState>>() = transition_rates;

        auto sim = mio::smm::Simulation(model, 0.0, 1.0);
        sim.advance(1.);
        //Save the number of transitions from region 0 to region 1
        transitions[n] += sim.get_model().populations[{mio::regions::Region(1), InfectionState::S}];
    }
    double mean_trans = std::accumulate(transitions.begin(), transitions.end(), 0.0) / num_runs;
    double rel_dev    = std::abs(expected_num_trans - mean_trans) / expected_num_trans;
    //The relative deviation should not exceed 15%
    EXPECT_LE(rel_dev, 0.15);
}
