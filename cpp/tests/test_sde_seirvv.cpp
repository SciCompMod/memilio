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
#include "sde_seirvv/model.h"
#include "sde_seirvv/simulation.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

const mio::sseirvv::Model<double>& sseirvv_testing_model()
{
    static mio::sseirvv::Model<double> model;
    model.step_size = 1.0 / 16;
    model.populations.array().setConstant(1);
    { // Set parameters s.t. coeffStoI is 1.
        model.parameters.set<mio::sseirvv::TimeExposedV1<double>>(1);
        model.parameters.set<mio::sseirvv::TimeExposedV2<double>>(1. / 4);
        model.parameters.set<mio::sseirvv::TimeInfectedV1<double>>(1);
        model.parameters.set<mio::sseirvv::TimeInfectedV2<double>>(1. / 4);
        model.parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV1<double>>(1);
        model.parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV2<double>>(1);
        model.parameters.get<mio::sseirvv::ContactPatterns<double>>().get_baseline()(0, 0) = 10;
    }
    return model;
}

TEST(TestSdeSeirvv, get_flows)
{
    // Make three get_flows computations with mocked rng.
    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<ScalarType>>>>>
        normal_dist_mock;

    // First two mock rng sets for test without clamping.
    // Third mock rng for test with clamping.
    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(27))
        .WillOnce(testing::Return(1.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(1.))
        .WillOnce(testing::Return(1.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(1.))
        .WillOnce(testing::Return(0.)) //End of first mock rng.
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(1.))
        .WillOnce(testing::Return(1.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(1.))
        .WillOnce(testing::Return(1.))
        .WillOnce(testing::Return(0.))
        .WillOnce(testing::Return(1.)) //End of second mock rng.
        .WillOnce(testing::Return(10.))
        .WillOnce(testing::Return(10.))
        .WillOnce(testing::Return(10.))
        .WillOnce(testing::Return(10.))
        .WillOnce(testing::Return(10.))
        .WillOnce(testing::Return(10.))
        .WillOnce(testing::Return(10.))
        .WillOnce(testing::Return(10.))
        .WillOnce(testing::Return(10.)); //End of third mock rng.

    // Non-constant y for a more meaningful test.
    Eigen::VectorXd y = Eigen::VectorXd(10);
    y << 1, 1, 1, 1, 1, 1, 1, 2, 2, 2;
    Eigen::VectorXd flows = Eigen::VectorXd::Constant(9, 1);

    // Results contain two parts : deterministic + stochastic.
    sseirvv_testing_model().get_flows(y, y, 0, flows);
    auto expected_result = Eigen::VectorXd(9);
    expected_result << 1 + sqrt(1) * sqrt(16), 3 + 0, 1 + 0, 4 + sqrt(4) * sqrt(16), 1 + sqrt(1) * sqrt(16), 4 + 0,
        3 + 0, 8 + sqrt(8) * sqrt(16), 8 + 0;
    EXPECT_EQ(flows, expected_result);

    sseirvv_testing_model().get_flows(y, y, 0, flows);
    expected_result = Eigen::VectorXd(9);
    expected_result << 1 + 0, 3 + sqrt(3) * sqrt(16), 1 + sqrt(1) * sqrt(16), 4 + 0, 1 + 0, 4 + sqrt(4) * sqrt(16),
        3 + sqrt(3) * sqrt(16), 8 + 0, 8 + sqrt(8) * sqrt(16);
    EXPECT_EQ(flows, expected_result);

    sseirvv_testing_model().get_flows(y, y, 0, flows);
    expected_result = Eigen::VectorXd(9);
    expected_result << 8, 8, 16, 16, 16, 16, 16, 32, 32;
    EXPECT_EQ(flows, expected_result);
}

TEST(TestSdeSeirvv, Simulation)
{
    // make a single integration step via a simulation
    // this should overwrite the model step size
    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<ScalarType>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(18))
        // 9 calls for each advance, as each call get_derivatives exactly once.
        .WillRepeatedly(testing::Return(.0));

    auto sim = mio::sseirvv::Simulation<double>(sseirvv_testing_model(), 0, 1);
    sim.advance(1);

    EXPECT_EQ(sim.get_model().step_size, 1.0); // Set by simulation.

    EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // Stores initial value and single step.

    auto expected_result = Eigen::VectorXd(10);
    expected_result << 0, 0.5, 1, 1, 0.5, 1, 2, 1, 1, 2;
    EXPECT_EQ(sim.get_result().get_last_value(), expected_result);

    sim.advance(1.5);
    EXPECT_EQ(sim.get_model().step_size, 0.5); // Set by simulation.
}

TEST(TestSdeSeirvv, FlowSimulation)
{
    // Make a single integration step via a flow simulation.
    // This should overwrite the model step size.
    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<ScalarType>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(18))
        // 3 calls for each advance, as each call get_derivatives exactly once.
        .WillRepeatedly(testing::Return(.5));

    auto sim = mio::sseirvv::FlowSimulation<double>(sseirvv_testing_model(), 0, 1);
    sim.advance(1);

    EXPECT_EQ(sim.get_model().step_size, 1.0); // Set by simulation.

    EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // Stores initial value and single step.

    auto expected_result = Eigen::VectorXd(10);
    expected_result << 0, 0.5, 1, 1, 0.5, 1, 2, 1, 1, 2;
    EXPECT_EQ(sim.get_result().get_last_value(), expected_result);

    auto expected_flows = Eigen::VectorXd(9);
    expected_flows << 0.5, 0.5, 1, 1, 1, 1, 1, 1, 1;
    EXPECT_EQ(sim.get_flows().get_last_value(), expected_flows);

    sim.advance(1.5);
    EXPECT_EQ(sim.get_model().step_size, 0.5); // Set by simulation.
}

TEST(TestSdeSeirvv, check_constraints_parameters)
{
    mio::sseirvv::Model<double>::ParameterSet parameters;
    parameters.set<mio::sseirvv::TimeInfectedV1<double>>(6);
    parameters.set<mio::sseirvv::TimeInfectedV2<double>>(6);
    parameters.set<mio::sseirvv::TimeExposedV1<double>>(6);
    parameters.set<mio::sseirvv::TimeExposedV2<double>>(6);
    parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV1<double>>(0.04);
    parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV2<double>>(0.04);
    parameters.get<mio::sseirvv::ContactPatterns<double>>().get_baseline()(0, 0) = 10;

    // model.check_constraints() combines the functions from population and parameters.
    // We only want to test the functions for the parameters defined in parameters.h.
    EXPECT_EQ(parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    parameters.set<mio::sseirvv::TimeInfectedV1<double>>(0);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::sseirvv::TimeInfectedV1<double>>(6);
    parameters.set<mio::sseirvv::TimeInfectedV2<double>>(0);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::sseirvv::TimeInfectedV2<double>>(6);
    parameters.set<mio::sseirvv::TimeExposedV1<double>>(0);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::sseirvv::TimeExposedV1<double>>(6);
    parameters.set<mio::sseirvv::TimeExposedV2<double>>(0);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::sseirvv::TimeExposedV2<double>>(6);
    parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV1<double>>(10.);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV1<double>>(0.04);
    parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV2<double>>(10.);
    EXPECT_EQ(parameters.check_constraints(), 1);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestSdeSeirvv, apply_constraints_parameters)
{
    const ScalarType tol_times = 1e-1;
    mio::sseirvv::Model<double>::ParameterSet parameters;
    parameters.set<mio::sseirvv::TimeInfectedV1<double>>(6);
    parameters.set<mio::sseirvv::TimeInfectedV2<double>>(6);
    parameters.set<mio::sseirvv::TimeExposedV1<double>>(6);
    parameters.set<mio::sseirvv::TimeExposedV2<double>>(6);
    parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV1<double>>(0.04);
    parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV2<double>>(0.04);
    parameters.get<mio::sseirvv::ContactPatterns<double>>().get_baseline()(0, 0) = 10;

    EXPECT_EQ(parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    parameters.set<mio::sseirvv::TimeInfectedV1<double>>(-2.5);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_EQ(parameters.get<mio::sseirvv::TimeInfectedV1<double>>(), tol_times);

    parameters.set<mio::sseirvv::TimeInfectedV2<double>>(-2.5);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_EQ(parameters.get<mio::sseirvv::TimeInfectedV2<double>>(), tol_times);

    parameters.set<mio::sseirvv::TimeExposedV1<double>>(-2.5);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_EQ(parameters.get<mio::sseirvv::TimeExposedV1<double>>(), tol_times);

    parameters.set<mio::sseirvv::TimeExposedV2<double>>(-2.5);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_EQ(parameters.get<mio::sseirvv::TimeExposedV2<double>>(), tol_times);

    parameters.set<mio::sseirvv::TimeExposedV2<double>>(-2.5);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_EQ(parameters.get<mio::sseirvv::TimeExposedV2<double>>(), tol_times);

    parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV1<double>>(10.);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_NEAR(parameters.get<mio::sseirvv::TransmissionProbabilityOnContactV1<double>>(), 0.0, 1e-14);

    parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV2<double>>(10.);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_NEAR(parameters.get<mio::sseirvv::TransmissionProbabilityOnContactV2<double>>(), 0.0, 1e-14);
    mio::set_log_level(mio::LogLevel::warn);
}
