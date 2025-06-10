/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Nils Wassmuth, Rene Schmieding, Martin J. Kuehn
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
#include "abm_helpers.h"
#include "sde_sir/model.h"
#include "sde_sir/simulation.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

const mio::ssir::Model<double>& ssir_testing_model()
{
    static mio::ssir::Model<double> model;
    model.step_size = 0.25;
    model.populations.array().setConstant(1);
    { // set parameters s.t. coeffStoI is 1
        model.parameters.set<mio::ssir::TimeInfected<double>>(1);
        model.parameters.set<mio::ssir::TransmissionProbabilityOnContact<double>>(1);
        model.parameters.get<mio::ssir::ContactPatterns<double>>().get_baseline()(0, 0) = 3;
    }
    return model;
}

TEST(TestSdeSir, get_flows)
{
    // make two get_flows computations with mocked rng
    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<double>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(4))
        .WillOnce(testing::Return(1.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(1.));

    const Eigen::Vector3d y = Eigen::Vector3d::Constant(1);
    Eigen::Vector2d flows   = Eigen::Vector2d::Constant(1);

    // results contain two parts : deterministic + stochastic

    ssir_testing_model().get_flows(y, y, 0, flows);
    auto expected_result = Eigen::Vector2d{1 + 2, 1 + 0};
    EXPECT_EQ(flows, expected_result);

    ssir_testing_model().get_flows(y, y, 0, flows);
    expected_result = Eigen::Vector2d{1 + 0, 1 + 2};
    EXPECT_EQ(flows, expected_result);
}

TEST(TestSdeSir, Simulation)
{
    // make a single integration step via a simulation
    // this should overwrite the model step size
    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<double>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(4))
        // 2 calls for each advance, as each call get_derivatives exactly once
        .WillRepeatedly(testing::Return(.5));

    auto sim = mio::ssir::Simulation(ssir_testing_model(), 0.0, 1.0);
    sim.advance(1.0);

    EXPECT_EQ(sim.get_model().step_size, 1.0); // set by simulation

    EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // stores initial value and single step

    auto expected_result = Eigen::Vector3d{0, 1, 2}; // flows each get capped to y[i]/step_size == 1/1
    EXPECT_EQ(sim.get_result().get_last_value(), expected_result);

    sim.advance(1.5);
    EXPECT_EQ(sim.get_model().step_size, 0.5); // set by simulation
}

TEST(TestSdeSir, FlowSimulation)
{
    // make a single integration step via a flow simulation
    // this should overwrite the model step size
    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<double>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(4))
        // 2 calls for each advance, as each call get_derivatives exactly once
        .WillRepeatedly(testing::Return(.5));

    auto sim = mio::ssir::FlowSimulation(ssir_testing_model(), 0.0, 1.0);
    sim.advance(1.0);

    EXPECT_EQ(sim.get_model().step_size, 1.0); // set by simulation

    EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // stores initial value and single step

    auto expected_result = Eigen::Vector3d{0, 1, 2}; // flows each get capped to y[i]/step_size == 1/1
    EXPECT_EQ(sim.get_result().get_last_value(), expected_result);

    auto expected_flows = Eigen::Vector2d{1, 1}; // flows each get capped to y[i]/step_size == 1/1
    EXPECT_EQ(sim.get_flows().get_last_value(), expected_flows);

    sim.advance(1.5);
    EXPECT_EQ(sim.get_model().step_size, 0.5); // set by simulation
}

TEST(TestSdeSir, check_constraints_parameters)
{
    mio::ssir::Model<double>::ParameterSet parameters;
    parameters.set<mio::ssir::TimeInfected<double>>(6);
    parameters.set<mio::ssir::TransmissionProbabilityOnContact<double>>(0.04);
    parameters.get<mio::ssir::ContactPatterns<double>>().get_baseline()(0, 0) = 10;

    // model.check_constraints() combines the functions from population and parameters.
    // We only want to test the functions for the parameters defined in parameters.h
    EXPECT_EQ(parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    parameters.set<mio::ssir::TimeInfected<double>>(0);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::ssir::TimeInfected<double>>(6);
    parameters.set<mio::ssir::TransmissionProbabilityOnContact<double>>(10.);
    EXPECT_EQ(parameters.check_constraints(), 1);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestSdeSir, apply_constraints_parameters)
{
    const double tol_times = 1e-1;
    mio::ssir::Model<double>::ParameterSet parameters;
    parameters.set<mio::ssir::TimeInfected<double>>(6);
    parameters.set<mio::ssir::TransmissionProbabilityOnContact<double>>(0.04);
    parameters.get<mio::ssir::ContactPatterns<double>>().get_baseline()(0, 0) = 10;

    EXPECT_EQ(parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    parameters.set<mio::ssir::TimeInfected<double>>(-2.5);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_EQ(parameters.get<mio::ssir::TimeInfected<double>>(), tol_times);

    parameters.set<mio::ssir::TransmissionProbabilityOnContact<double>>(10.);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_NEAR(parameters.get<mio::ssir::TransmissionProbabilityOnContact<double>>(), 0.0, 1e-14);
    mio::set_log_level(mio::LogLevel::warn);
}
