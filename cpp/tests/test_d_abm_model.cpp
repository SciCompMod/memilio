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

#include "d_abm/parameters.h"
#include "d_abm/quad_well.h"
#include "d_abm/model.h"
#include "d_abm/simulation.h"
#include "memilio/utils/random_number_generator.h"
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
    auto well = well_index(Eigen::Vector2d{-1, 1});
    EXPECT_EQ(well, size_t(0));
}

TEST(TestQuadWell, adopt)
{
    QuadWellModel<InfectionState>::Agent agent{Eigen::Vector2d{-1, 1}, InfectionState::S};
    QuadWellModel<InfectionState> qw({agent}, {});
    EXPECT_EQ(agent.status, InfectionState::S);
    qw.adopt(agent, InfectionState::E);
    EXPECT_EQ(agent.status, InfectionState::E);
}

TEST(TestQuadWell, adoptionRate)
{
    QuadWellModel<InfectionState>::Agent a1{Eigen::Vector2d{-1, 1}, InfectionState::S};
    QuadWellModel<InfectionState>::Agent a2{Eigen::Vector2d{-1.2, 1}, InfectionState::I};
    QuadWellModel<InfectionState> qw({a1, a2}, {{InfectionState::S,
                                                 InfectionState::E,
                                                 mio::dabm::Region(0),
                                                 0.1,
                                                 {InfectionState::C, InfectionState::I},
                                                 {1, 0.5}}});
    EXPECT_EQ(qw.adoption_rate(a1, InfectionState::E), 0.025);
}

TEST(TestQuadWell, move)
{
    QuadWellModel<InfectionState>::Agent a1{Eigen::Vector2d{-1.2, 1}, InfectionState::S};
    QuadWellModel<InfectionState>::Agent a2{Eigen::Vector2d{-1.2, 1}, InfectionState::I};
    QuadWellModel<InfectionState> qw({a1, a2}, {}, 0.1, 0., {InfectionState::I});
    qw.move(0, 0.1, a1);
    qw.move(0, 0.1, a2);
    EXPECT_EQ(a1.position[0], -0.9888);
    EXPECT_EQ(a1.position[1], 1);
    EXPECT_EQ(a2.position[0], -1.2);
    EXPECT_EQ(a2.position[1], 1);
}

TEST(TestQuadWell, setNonMovingRegions)
{
    QuadWellModel<InfectionState>::Agent a1{Eigen::Vector2d{-1.2, 1}, InfectionState::S};
    QuadWellModel<InfectionState> qw({a1}, {});
    qw.set_non_moving_regions({0});
    qw.move(0, 0.1, a1);
    EXPECT_EQ(a1.position[0], -1.2);
    EXPECT_EQ(a1.position[1], 1);
}

TEST(TestDABMSimulation, advance)
{
    using testing::Return;
    using Model = mio::dabm::Model<QuadWellModel<InfectionState>>;
    QuadWellModel<InfectionState>::Agent a1{Eigen::Vector2d{-1, 1}, InfectionState::S};
    QuadWellModel<InfectionState>::Agent a2{Eigen::Vector2d{-1, 1}, InfectionState::R};
    QuadWellModel<InfectionState>::Agent a3{Eigen::Vector2d{-1, 1}, InfectionState::I};
    std::vector<mio::dabm::AdoptionRate<InfectionState>> adoption_rates;
    for (size_t region = 0; region < 4; ++region) {
        adoption_rates.push_back({InfectionState::S,
                                  InfectionState::E,
                                  mio::dabm::Region(region),
                                  0.1,
                                  {InfectionState::C, InfectionState::I},
                                  {1, 0.5}});
        adoption_rates.push_back({InfectionState::E, InfectionState::C, mio::dabm::Region(region), 1.0 / 5., {}, {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::R, mio::dabm::Region(region), 0.2 / 3., {}, {}});
        adoption_rates.push_back({InfectionState::C, InfectionState::I, mio::dabm::Region(region), 0.8 / 3., {}, {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::R, mio::dabm::Region(region), 0.99 / 5., {}, {}});
        adoption_rates.push_back({InfectionState::I, InfectionState::D, mio::dabm::Region(region), 0.01 / 5., {}, {}});
    }
    Model model({a1, a2, a3}, adoption_rates, 0.4, 0.0, {InfectionState::D});
    auto sim = mio::Simulation<double, Model>(model, 0.0, 0.1);
    //Setup so first adoption event will be in second time step
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::AtLeast(1))
        .WillOnce(Return(0.0226))
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
