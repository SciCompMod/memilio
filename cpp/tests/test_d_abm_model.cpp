/* 
* Copyright (C) 2020-2026 MEmilio
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

#include "d_abm/parameters.h"
#include "d_abm/quad_well.h"
#include "d_abm/single_well.h"
#include "d_abm/model.h"
#include "d_abm/simulation.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/geography/regions.h"
#include "abm_helpers.h"
#include <cstddef>
#include <gtest/gtest.h>
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

TEST(TestQuadWell, well)
{
    //Test whether the well indices are calculated correctly
    EXPECT_EQ(well_index(Eigen::Vector2d{-1, 1}), size_t(0));
    EXPECT_EQ(well_index(Eigen::Vector2d{1, 1}), size_t(1));
    EXPECT_EQ(well_index(Eigen::Vector2d{-1, -1}), size_t(2));
    EXPECT_EQ(well_index(Eigen::Vector2d{1, -1}), size_t(3));
}

TEST(TestQuadWell, adopt)
{
    //Test adopt function i.e. if agent adopts the given status
    QuadWell<InfectionState>::Agent agent{Eigen::Vector2d{-1, 1}, InfectionState::S};
    QuadWell<InfectionState> qw({agent}, {});
    EXPECT_EQ(agent.status, InfectionState::S);
    qw.adopt(agent, InfectionState::E);
    EXPECT_EQ(agent.status, InfectionState::E);
}

TEST(TestSingleWell, adopt)
{
    //Test adopt function i.e. if agent adopts the given status
    SingleWell<InfectionState>::Agent agent{Eigen::Vector2d{-1, 1}, InfectionState::S};
    SingleWell<InfectionState> sw({agent}, {});
    EXPECT_EQ(agent.status, InfectionState::S);
    sw.adopt(agent, InfectionState::E);
    EXPECT_EQ(agent.status, InfectionState::E);
}

TEST(TestQuadWell, adoptionRate)
{
    //Test whether adoption rates are calculated correctly.
    //First-order adoption rates are given by:
    // rate.factor * N(rate.from, rate.region)
    // with N(from, region) the number of agents in Region "region" having infection state "from"
    //Second-order adoption rates are given by:
    // rate.factor * N(rate.from, rate.region)/total_pop * sum (over all rate.influences) influence.factor * N_contact(influence.status, rate.region)
    // with N_contact(status, region) the number of agents in Region "region" having infection state "status" that are within the contact radius

    //Agents in region 0
    QuadWell<InfectionState>::Agent a1{Eigen::Vector2d{-1, 1}, InfectionState::S};
    QuadWell<InfectionState>::Agent a2{Eigen::Vector2d{-1.2, 1}, InfectionState::I};
    //Agent in region 1
    QuadWell<InfectionState>::Agent a3{Eigen::Vector2d{1, 1}, InfectionState::I};
    //Agent in region 0
    QuadWell<InfectionState>::Agent a4{Eigen::Vector2d{-1.1, 1}, InfectionState::I};
    //Initialize model without fourth agent
    QuadWell<InfectionState> qw({a1, a2, a3}, {{InfectionState::S,
                                                InfectionState::E,
                                                mio::regions::Region(0),
                                                0.1,
                                                {{InfectionState::C, 1}, {InfectionState::I, 0.5}}}});
    //Initialize model with all agents
    QuadWell<InfectionState> qw1({a1, a2, a3, a4}, {{InfectionState::S,
                                                     InfectionState::E,
                                                     mio::regions::Region(0),
                                                     0.1,
                                                     {{InfectionState::C, 1}, {InfectionState::I, 0.5}}}});
    //a1 only has contact to a2 as a3 is in another region
    EXPECT_EQ(qw.adoption_rate(a1, InfectionState::E), 0.025);
    //There is no rate from I to E, hence rate for a2 is 0
    EXPECT_EQ(qw.adoption_rate(a2, InfectionState::E), 0.0);
    //a1 now has contact to a2 and a4 (both status I) which increases the rate
    EXPECT_EQ(qw1.adoption_rate(a1, InfectionState::E), 1. / 30.);
}

TEST(TestSingleWell, adoptionRate)
{
    //Test whether adoption rates are calculated correctly.
    //First-order adoption rates are given by:
    // rate.factor * N(rate.from)
    // with N(from) the number of agents having infection state "from"
    // Second-order adoption rates are given by:
    // rate.factor * N(rate.from)/total_pop * sum (over all rate.influences) influence.factor * N_contact(influence.status)
    // with N_contact(status) the number of agents having infection state "status" that are within the contact radius

    //Agents
    SingleWell<InfectionState>::Agent a1{Eigen::Vector2d{-1, 1}, InfectionState::S};
    SingleWell<InfectionState>::Agent a2{Eigen::Vector2d{-1.2, 1}, InfectionState::I};
    SingleWell<InfectionState>::Agent a3{Eigen::Vector2d{-1.1, 1}, InfectionState::I};
    //Initialize model without third agent
    SingleWell<InfectionState> sw({a1, a2},
                                  {{InfectionState::S,
                                    InfectionState::E,
                                    mio::regions::Region(0),
                                    0.1,
                                    {{InfectionState::C, 1}, {InfectionState::I, 0.5}}},
                                   {InfectionState::I, InfectionState::R, mio::regions::Region(0), 0.15, {}}});
    //Initialize model with all agents
    SingleWell<InfectionState> sw1({a1, a2, a3}, {{InfectionState::S,
                                                   InfectionState::E,
                                                   mio::regions::Region(0),
                                                   0.1,
                                                   {{InfectionState::C, 1}, {InfectionState::I, 0.5}}}});
    //a1 has contact to a2
    EXPECT_EQ(sw.adoption_rate(a1, InfectionState::E), 0.025);
    //a2 can recover
    EXPECT_EQ(sw.adoption_rate(a2, InfectionState::R), 0.15);
    //There is no rate from I to E, hence rate for a2 is 0
    EXPECT_EQ(sw.adoption_rate(a2, InfectionState::E), 0.0);
    //a1 now has contact to a2 and a3 (both status I) which increases the rate
    EXPECT_EQ(sw1.adoption_rate(a1, InfectionState::E), 1. / 30.);
}

TEST(TestQuadWell, move)
{
    //Test evaluation of diffusion process with euler-maruyama integration
    QuadWell<InfectionState>::Agent a1{Eigen::Vector2d{-1.2, 1}, InfectionState::S};
    QuadWell<InfectionState>::Agent a2{Eigen::Vector2d{-1.2, 1}, InfectionState::I};
    //Sigma is set to 0, thus movement is only given by the function gradient
    //Set I as non-moving state
    QuadWell<InfectionState> qw({a1, a2}, {}, 0.1, 0., {InfectionState::I});
    qw.move(0, 0.1, a1);
    qw.move(0, 0.1, a2);
    //a1 moves in direction of the function gradient
    EXPECT_NEAR(a1.position[0], -0.9888, 1e-12);
    EXPECT_NEAR(a1.position[1], 1, 1e-12);
    //As a2 has infection state I, it does not move
    EXPECT_EQ(a2.position[0], -1.2);
    EXPECT_EQ(a2.position[1], 1);
}

TEST(TestSingleWell, move)
{
    //Test evaluation of diffusion process with euler-maruyama integration
    SingleWell<InfectionState>::Agent a1{Eigen::Vector2d{-1.2, 1}, InfectionState::S};
    SingleWell<InfectionState>::Agent a2{Eigen::Vector2d{-1.2, 1}, InfectionState::I};
    //Sigma is set to 0, thus movement is only given by the function gradient
    //Set I as non-moving state
    SingleWell<InfectionState> sw({a1, a2}, {}, 0.1, 0., {InfectionState::I});
    sw.move(0, 0.1, a1);
    sw.move(0, 0.1, a2);
    //a1 moves in direction of the function gradient
    EXPECT_NEAR(a1.position[0], -0.8544, 1e-12);
    EXPECT_NEAR(a1.position[1], 0.8, 1e-12);
    //As a2 has infection state I, it does not move
    EXPECT_EQ(a2.position[0], -1.2);
    EXPECT_EQ(a2.position[1], 1);
}

TEST(TestSingleWell, time_point)
{
    //Test evaluation time_point function of singlewell potential
    SingleWell<InfectionState>::Agent a1{Eigen::Vector2d{-1.2, 1}, InfectionState::S};
    SingleWell<InfectionState>::Agent a2{Eigen::Vector2d{-1.2, 1}, InfectionState::I};
    SingleWell<InfectionState>::Agent a3{Eigen::Vector2d{1, 1}, InfectionState::I};
    SingleWell<InfectionState> sw({a1, a2, a3}, {}, 0.1, 0., {});
    auto vec = sw.time_point();
    //Agents are aggregated by their compartment
    EXPECT_EQ(vec.size(), static_cast<int>(InfectionState::Count));
    EXPECT_NEAR(vec[static_cast<size_t>(InfectionState::S)], 1, 1e-14);
    EXPECT_NEAR(vec[static_cast<size_t>(InfectionState::E)], 0, 1e-14);
    EXPECT_NEAR(vec[static_cast<size_t>(InfectionState::I)], 2, 1e-14);
}

TEST(TestQuadWell, change_well)
{
    //Test spatial transition counting
    QuadWell<InfectionState>::Agent a{Eigen::Vector2d{-0.1, 1}, InfectionState::S};
    //Sigma is set to 2
    QuadWell<InfectionState> qw({a}, {}, 0.1, 2.);
    //mock value for Brownian motion
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::NormalDistribution<double>>>> mock_normal_dist;
    EXPECT_CALL(mock_normal_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.5)) // noise in x-direction
        .WillRepeatedly(testing::Return(0.0));
    qw.move(0, 0.1, a);
    //a changes from well 0 to well 1
    EXPECT_EQ(well_index(a.position), size_t(1));
    EXPECT_EQ(qw.number_transitions()[static_cast<int>(InfectionState::S)](0, 1), 1);
}

TEST(TestQuadWell, setNonMovingRegions)
{
    //Test non-moving regions
    //a1 is in region 0
    QuadWell<InfectionState>::Agent a1{Eigen::Vector2d{-1.2, 1}, InfectionState::S};
    QuadWell<InfectionState> qw({a1}, {});
    //Set region 0 as non-moving regions
    qw.set_non_moving_regions({0});
    qw.move(0, 0.1, a1);
    EXPECT_EQ(a1.position[0], -1.2);
    EXPECT_EQ(a1.position[1], 1);
}

TEST(TestDABMSimulation, advance)
{
    //Test whether temporal Gillespie algorithm calculates events in the correct order
    using testing::Return;
    using Model = mio::dabm::Model<QuadWell<InfectionState>>;
    QuadWell<InfectionState>::Agent a1{Eigen::Vector2d{-1, 1}, InfectionState::S};
    QuadWell<InfectionState>::Agent a2{Eigen::Vector2d{-1, 1}, InfectionState::R};
    QuadWell<InfectionState>::Agent a3{Eigen::Vector2d{-1, 1}, InfectionState::I};
    std::vector<mio::AdoptionRate<double, InfectionState>> adoption_rates;
    //Add adoption rates for every region
    for (size_t region = 0; region < 4; ++region) {
        adoption_rates.push_back({InfectionState::S,
                                  InfectionState::E,
                                  mio::regions::Region(region),
                                  0.1,
                                  {{InfectionState::C, 1}, {InfectionState::I, 0.5}}});
        adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::regions::Region(region), 1.0 / 5., {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::R, mio::regions::Region(region), 0.2 / 3., {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::I, mio::regions::Region(region), 0.8 / 3., {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::R, mio::regions::Region(region), 0.99 / 5., {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::D, mio::regions::Region(region), 0.01 / 5., {}});
    }
    Model model({a1, a2, a3}, adoption_rates, 0.4, 0.0, {InfectionState::D});
    auto sim = mio::dabm::Simulation(model, 0.0, 0.1);
    //Setup such that first adoption event will be in second time step
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::AtLeast(1))
        .WillOnce(Return(0.0226)) // Waiting time for first adoption rate
        .WillRepeatedly(testing::Return(1.0));
    //Setup so first adoption event will be the transition of a3 from I to R
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(size_t(1)));
    sim.advance(30.);
    //Initial positions and infection state of all agents
    EXPECT_EQ(sim.get_result().get_value(0)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(0)[static_cast<size_t>(InfectionState::I)], 1);
    EXPECT_EQ(sim.get_result().get_value(0)[static_cast<size_t>(InfectionState::R)], 1);
    //In the first time step no adoption event happens
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::I)], 1);
    EXPECT_EQ(sim.get_result().get_value(1)[static_cast<size_t>(InfectionState::R)], 1);
    //a3 transitions from I to R in the second time step
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::S)], 1);
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::I)], 0);
    EXPECT_EQ(sim.get_result().get_value(2)[static_cast<size_t>(InfectionState::R)], 2);
    //check whether simulation advances until the end
    EXPECT_EQ(sim.get_result().get_last_time(), 30.);
}
