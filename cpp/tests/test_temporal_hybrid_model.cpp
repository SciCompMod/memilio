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

#include "d_abm/model.h"
#include "d_abm/simulation.h"
#include "d_abm/single_well.h"
#include "hybrid/infection_state.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/geography/regions.h"
#include "memilio/utils/random_number_generator.h"
#include "models/hybrid/conversion_functions.cpp"
#include "models/hybrid/temporal_hybrid_model.h"
#include "ode_secir/infection_state.h"
#include "smm/model.h"
#include "smm/simulation.h"
#include "abm_helpers.h"
#include <gtest/gtest.h>
#include <vector>

struct MockModel1 {
    double a1;
    double a2;
    std::vector<double> result;

    void update_result()
    {
        result = {a1, a2};
    }

    void advance(double /*tmax*/)
    {
        a1 += 0.3;
        a2 *= 2.;
        update_result();
    }
};

struct MockModel2 {
    double pop;
    double result;

    void update_result()
    {
        result = pop;
    }

    void advance(double /*tmax*/)
    {
        if (pop > 1) {
            pop--;
        }
        else {
            pop = 1.5;
        }
        update_result();
    }
};

template <>
void mio::hybrid::convert_model(const MockModel1& model1, MockModel2& model2)
{
    model2.pop = model1.a1 + model1.a2;
    model2.update_result();
}

template <>
void mio::hybrid::convert_model(const MockModel2& model2, MockModel1& model1)
{
    model1.a1 = model2.pop / 2.;
    model1.a2 = model2.pop / 2.;
    model1.update_result();
}

std::vector<double> result_fct_m1(const MockModel1& model1, double /*t*/)
{
    return model1.result;
}

double result_fct_m2(const MockModel2& model2, double /*t*/)
{
    return model2.result;
}

/**
 * @brief Test advance function of TemporalHybridSimulation class.
 * The advance function of this class should advance one of two models depending on a given condition. 
 * If the condition is fulfilled, the currently used model used be converted to the target model using the corresponding convert_model function.
 */
TEST(TestTemporalHybrid, test_advance)
{
    //Create two models that should be combined in a temporl-hybrid model
    MockModel1 m1{0.1, 0.4, {}};
    m1.update_result();
    MockModel2 m2;

    mio::hybrid::TemporalHybridSimulation<MockModel1, MockModel2, std::vector<double>, double> hybrid_sim(
        m1, m2, &result_fct_m1, &result_fct_m2, true, 0., 0.5);

    //If the total population is bigger than 1, model2 should be used and otherwise model1 should be used
    const auto condition = [](const std::vector<double>& result_m1, const double result_m2, bool model1_used) {
        if (model1_used) {
            if (std::accumulate(result_m1.begin(), result_m1.end(), 0.0) > 1) {
                return true;
            }
        }
        else {
            if (result_m2 < 1) {
                return true;
            }
        }
        return false;
    };

    hybrid_sim.advance(0.5, condition);

    EXPECT_EQ(hybrid_sim.using_model1(), true);
    EXPECT_NEAR(hybrid_sim.get_model1().a1, 0.4, 1e-10);
    EXPECT_NEAR(hybrid_sim.get_model1().a2, 0.8, 1e-10);

    hybrid_sim.advance(1., condition);

    EXPECT_EQ(hybrid_sim.using_model1(), false);
    EXPECT_NEAR(hybrid_sim.get_model2().pop, 0.2, 1e-10);

    hybrid_sim.advance(1.5, condition);

    EXPECT_EQ(hybrid_sim.using_model1(), true);
    EXPECT_NEAR(hybrid_sim.get_model1().a1, 0.4, 1e-10);
    EXPECT_NEAR(hybrid_sim.get_model1().a2, 0.2, 1e-10);
}

/**
 * @brief Test conversion from dABM to SMM and vice versa.
 */
TEST(TestTemporalHybrid, test_conversion_dabm_smm)
{
    using Model1 = mio::dabm::Model<SingleWell<mio::hybrid::InfectionState>>;
    using Model2 = mio::smm::Model<1, mio::hybrid::InfectionState>;

    //Initialize agents for dabm
    SingleWell<mio::hybrid::InfectionState>::Agent a1{Eigen::Vector2d{-0.5, 0},
                                                      mio::hybrid::InfectionState::Susceptible};
    SingleWell<mio::hybrid::InfectionState>::Agent a2{Eigen::Vector2d{0.5, 0},
                                                      mio::hybrid::InfectionState::Susceptible};
    SingleWell<mio::hybrid::InfectionState>::Agent a3{Eigen::Vector2d{0.5, 0.5},
                                                      mio::hybrid::InfectionState::InfectedSymptoms};

    Model1 model1({a1, a2, a3}, {});
    Model2 model2;
    model2.parameters.get<mio::smm::AdoptionRates<mio::hybrid::InfectionState>>().push_back(
        {mio::hybrid::InfectionState::Susceptible,
         mio::hybrid::InfectionState::Exposed,
         mio::regions::Region(0),
         0.1,
         {{mio::hybrid::InfectionState::InfectedNoSymptoms, 1}, {mio::hybrid::InfectionState::InfectedSymptoms, 0.5}}});

    //Parameters for simulation
    double t0 = 0;
    double dt = 0.1;

    auto sim_dabm = mio::dabm::Simulation(model1, t0, dt);
    auto sim_smm  = mio::smm::Simulation(model2, t0 - 1, dt);

    //Convert dabm simulation to smm simulation
    mio::hybrid::convert_model(sim_dabm, sim_smm);

    EXPECT_EQ(sim_smm.get_result().get_last_time(), t0);
    EXPECT_NEAR(sim_smm.get_result().get_last_value()[(int)mio::hybrid::InfectionState::Susceptible], 2, 1e-10);
    EXPECT_NEAR(sim_smm.get_result().get_last_value()[(int)mio::hybrid::InfectionState::InfectedSymptoms], 1, 1e-10);
    auto pop_S = sim_smm.get_model().populations[{mio::regions::Region(0), mio::hybrid::InfectionState::Susceptible}];
    auto pop_ISy =
        sim_smm.get_model().populations[{mio::regions::Region(0), mio::hybrid::InfectionState::InfectedSymptoms}];
    EXPECT_NEAR(pop_S, 2, 1e-10);
    EXPECT_NEAR(pop_ISy, 1, 1e-10);

    //Delete dabm population
    sim_dabm.get_model().populations.clear();

    EXPECT_EQ(sim_dabm.get_model().populations.size(), 0);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    //Distribution to sample agents' position
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(6))
        .WillOnce(testing::Return(-0.5)) // x-value agent1
        .WillOnce(testing::Return(0)) // y-value agent1
        .WillOnce(testing::Return(0.5)) // x-value agent2
        .WillOnce(testing::Return(0)) // y-value agent2
        .WillOnce(testing::Return(0.5)) // x-value agent3
        .WillOnce(testing::Return(0.5)) // y-value agent3
        .WillRepeatedly(testing::Return(1.0));

    //Distribution to sample agents' infection state
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke)
        .Times(testing::AtLeast(3))
        .WillOnce(testing::Return(0)) // agent1
        .WillOnce(testing::Return(0)) // agent2
        .WillOnce(testing::Return(3)) // agent3
        .WillRepeatedly(testing::Return(1));

    //Convert smm simulation to dabm simulation
    mio::hybrid::convert_model(sim_smm, sim_dabm);

    EXPECT_EQ(sim_dabm.get_model().populations.size(), 3);
    //agent1
    EXPECT_EQ(sim_dabm.get_model().populations[0].position[0], -0.5);
    EXPECT_EQ(sim_dabm.get_model().populations[0].position[1], 0);
    EXPECT_EQ(sim_dabm.get_model().populations[0].status, mio::hybrid::InfectionState::Susceptible);
    //agent2
    EXPECT_EQ(sim_dabm.get_model().populations[1].position[0], 0.5);
    EXPECT_EQ(sim_dabm.get_model().populations[1].position[1], 0);
    EXPECT_EQ(sim_dabm.get_model().populations[1].status, mio::hybrid::InfectionState::Susceptible);
    //agent3
    EXPECT_EQ(sim_dabm.get_model().populations[2].position[0], 0.5);
    EXPECT_EQ(sim_dabm.get_model().populations[2].position[1], 0.5);
    EXPECT_EQ(sim_dabm.get_model().populations[2].status, mio::hybrid::InfectionState::InfectedSymptoms);

    //Test if conversion also works if agents should just be overwritten
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist1;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>>
        mock_discrete_dist1;

    //Distribution to sample agents' position
    EXPECT_CALL(mock_uniform_dist1.get_mock(), invoke)
        .Times(testing::AtLeast(6))
        .WillOnce(testing::Return(-0.1)) // x-value agent1
        .WillOnce(testing::Return(0.1)) // y-value agent1
        .WillOnce(testing::Return(0.1)) // x-value agent2
        .WillOnce(testing::Return(-0.1)) // y-value agent2
        .WillOnce(testing::Return(0.1)) // x-value agent3
        .WillOnce(testing::Return(0.1)) // y-value agent3
        .WillRepeatedly(testing::Return(1.0));

    //Distribution to sample agents' infection state
    EXPECT_CALL(mock_discrete_dist1.get_mock(), invoke)
        .Times(testing::AtLeast(3))
        .WillOnce(testing::Return(2)) // agent1
        .WillOnce(testing::Return(2)) // agent2
        .WillOnce(testing::Return(0)) // agent3
        .WillRepeatedly(testing::Return(1));

    //Convert smm simulation to dabm simulation
    mio::hybrid::convert_model(sim_smm, sim_dabm);

    EXPECT_EQ(sim_dabm.get_model().populations.size(), 3);
    //agent1
    EXPECT_EQ(sim_dabm.get_model().populations[0].position[0], -0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[0].position[1], 0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[0].status, mio::hybrid::InfectionState::InfectedNoSymptoms);
    //agent2
    EXPECT_EQ(sim_dabm.get_model().populations[1].position[0], 0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[1].position[1], -0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[1].status, mio::hybrid::InfectionState::InfectedNoSymptoms);
    //agent3
    EXPECT_EQ(sim_dabm.get_model().populations[2].position[0], 0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[2].position[1], 0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[2].status, mio::hybrid::InfectionState::Susceptible);
}

/**
 * @brief Test conversion from dABM to ode-secir and vice versa.
 */
TEST(TestTemporalHybrid, test_conversion_dabm_osecir)
{
    using Model1 = mio::dabm::Model<SingleWell<mio::hybrid::InfectionState>>;
    using Model2 = mio::osecir::Model<double>;

    //Initialize agents for dabm
    SingleWell<mio::hybrid::InfectionState>::Agent a1{Eigen::Vector2d{-0.5, 0},
                                                      mio::hybrid::InfectionState::Susceptible};
    SingleWell<mio::hybrid::InfectionState>::Agent a2{Eigen::Vector2d{0.5, 0},
                                                      mio::hybrid::InfectionState::Susceptible};
    SingleWell<mio::hybrid::InfectionState>::Agent a3{Eigen::Vector2d{0.5, 0.5},
                                                      mio::hybrid::InfectionState::InfectedSymptoms};

    Model1 model1({a1, a2, a3}, {});
    Model2 model2(2);

    //Parameters for simulation
    double t0 = 0;
    double dt = 0.1;

    auto sim_dabm   = mio::dabm::Simulation(model1, t0, dt);
    auto sim_osecir = mio::Simulation(model2, t0 - 1, dt);

    //Convert dabm simulation to osecir simulation
    mio::hybrid::convert_model(sim_dabm, sim_osecir);

    EXPECT_EQ(sim_osecir.get_result().get_last_time(), t0);
    EXPECT_NEAR(sim_osecir.get_result().get_last_value()[sim_osecir.get_model().populations.get_flat_index(
                    {mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible})],
                1, 1e-10);
    EXPECT_NEAR(sim_osecir.get_result().get_last_value()[sim_osecir.get_model().populations.get_flat_index(
                    {mio::AgeGroup(1), mio::osecir::InfectionState::Susceptible})],
                1, 1e-10);
    EXPECT_NEAR(sim_osecir.get_result().get_last_value()[sim_osecir.get_model().populations.get_flat_index(
                    {mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms})],
                0.5, 1e-10);
    EXPECT_NEAR(sim_osecir.get_result().get_last_value()[sim_osecir.get_model().populations.get_flat_index(
                    {mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSymptoms})],
                0.5, 1e-10);
    auto pop_S = sim_osecir.get_model().populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] +
                 sim_osecir.get_model().populations[{mio::AgeGroup(1), mio::osecir::InfectionState::Susceptible}];
    auto pop_ISy =
        sim_osecir.get_model().populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] +
        sim_osecir.get_model().populations[{mio::AgeGroup(1), mio::osecir::InfectionState::InfectedSymptoms}];
    EXPECT_NEAR(pop_S, 2, 1e-10);
    EXPECT_NEAR(pop_ISy, 1, 1e-10);

    //Delete dabm population
    sim_dabm.get_model().populations.clear();

    EXPECT_EQ(sim_dabm.get_model().populations.size(), 0);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;

    //Distribution to sample agents' position
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(6))
        .WillOnce(testing::Return(-0.5)) // x-value agent1
        .WillOnce(testing::Return(0)) // y-value agent1
        .WillOnce(testing::Return(0.5)) // x-value agent2
        .WillOnce(testing::Return(0)) // y-value agent2
        .WillOnce(testing::Return(0.5)) // x-value agent3
        .WillOnce(testing::Return(0.5)) // y-value agent3
        .WillRepeatedly(testing::Return(1.0));

    //Distribution to sample agents' infection state
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke)
        .Times(testing::AtLeast(3))
        .WillOnce(testing::Return(0)) // agent1
        .WillOnce(testing::Return(0)) // agent2
        .WillOnce(testing::Return(3)) // agent3
        .WillRepeatedly(testing::Return(1));

    //Convert ode-secir simulation to dabm simulation
    mio::hybrid::convert_model(sim_osecir, sim_dabm);

    EXPECT_EQ(sim_dabm.get_model().populations.size(), 3);
    //agent1
    EXPECT_EQ(sim_dabm.get_model().populations[0].position[0], -0.5);
    EXPECT_EQ(sim_dabm.get_model().populations[0].position[1], 0);
    EXPECT_EQ(sim_dabm.get_model().populations[0].status, mio::hybrid::InfectionState::Susceptible);
    //agent2
    EXPECT_EQ(sim_dabm.get_model().populations[1].position[0], 0.5);
    EXPECT_EQ(sim_dabm.get_model().populations[1].position[1], 0);
    EXPECT_EQ(sim_dabm.get_model().populations[1].status, mio::hybrid::InfectionState::Susceptible);
    //agent3
    EXPECT_EQ(sim_dabm.get_model().populations[2].position[0], 0.5);
    EXPECT_EQ(sim_dabm.get_model().populations[2].position[1], 0.5);
    EXPECT_EQ(sim_dabm.get_model().populations[2].status, mio::hybrid::InfectionState::InfectedSymptoms);

    //Test if conversion also works if agents should just be overwritten
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist1;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>>
        mock_discrete_dist1;

    //Distribution to sample agents' position
    EXPECT_CALL(mock_uniform_dist1.get_mock(), invoke)
        .Times(testing::AtLeast(6))
        .WillOnce(testing::Return(-0.1)) // x-value agent1
        .WillOnce(testing::Return(0.1)) // y-value agent1
        .WillOnce(testing::Return(0.1)) // x-value agent2
        .WillOnce(testing::Return(-0.1)) // y-value agent2
        .WillOnce(testing::Return(0.1)) // x-value agent3
        .WillOnce(testing::Return(0.1)) // y-value agent3
        .WillRepeatedly(testing::Return(1.0));

    //Distribution to sample agents' infection state
    EXPECT_CALL(mock_discrete_dist1.get_mock(), invoke)
        .Times(testing::AtLeast(3))
        .WillOnce(testing::Return(2)) // agent1
        .WillOnce(testing::Return(2)) // agent2
        .WillOnce(testing::Return(0)) // agent3
        .WillRepeatedly(testing::Return(1));

    //Convert ode-secir simulation to dabm simulation
    mio::hybrid::convert_model(sim_osecir, sim_dabm);

    EXPECT_EQ(sim_dabm.get_model().populations.size(), 3);
    //agent1
    EXPECT_EQ(sim_dabm.get_model().populations[0].position[0], -0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[0].position[1], 0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[0].status, mio::hybrid::InfectionState::InfectedNoSymptoms);
    //agent2
    EXPECT_EQ(sim_dabm.get_model().populations[1].position[0], 0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[1].position[1], -0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[1].status, mio::hybrid::InfectionState::InfectedNoSymptoms);
    //agent3
    EXPECT_EQ(sim_dabm.get_model().populations[2].position[0], 0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[2].position[1], 0.1);
    EXPECT_EQ(sim_dabm.get_model().populations[2].status, mio::hybrid::InfectionState::Susceptible);
}
