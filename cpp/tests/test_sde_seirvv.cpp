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
#include "memilio/compartments/stochastic_simulation.h"
#include "memilio/config.h"
#include "sde_seirvv/model.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

const mio::sseirvv::Model<double>& sseirvv_testing_model()
{
    static mio::sseirvv::Model<double> model;
    model.populations.array().setConstant(1);
    {
        model.parameters.set<mio::sseirvv::TimeExposedV1<double>>(2);
        model.parameters.set<mio::sseirvv::TimeExposedV2<double>>(4);
        model.parameters.set<mio::sseirvv::TimeInfectedV1<double>>(4);
        model.parameters.set<mio::sseirvv::TimeInfectedV2<double>>(2);
        model.parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV1<double>>(0.5);
        model.parameters.set<mio::sseirvv::TransmissionProbabilityOnContactV2<double>>(0.25);
        model.parameters.get<mio::sseirvv::ContactPatterns<double>>().get_baseline()(0, 0) = 10;
    }
    return model;
}

TEST(TestSdeSeirvv, Model)
{
    // check get_flows and get_noise
    const Eigen::VectorXd y = Eigen::VectorXd::Constant(10, 1);
    Eigen::VectorXd flows   = Eigen::VectorXd::Constant(9, 1);
    Eigen::VectorXd noise   = Eigen::VectorXd::Constant(10, 1);

    sseirvv_testing_model().get_flows(y, y, 0, flows);
    auto expected_flows = Eigen::Vector<double, 9>{0.5, 0.5, 0.5, 0.25, 0.25, 0.5, 0.5, 0.25, 0.5};
    EXPECT_EQ(flows, expected_flows);

    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<double>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(9)) // one call for each flow
        .WillRepeatedly(testing::Return(0.5));

    sseirvv_testing_model().get_noise(y, y, 0, noise);
    Eigen::VectorXd expected_noise(10);
    sseirvv_testing_model().get_derivatives(expected_flows.array().sqrt() * 0.5, expected_noise);
    EXPECT_EQ(noise, expected_noise);
}

TEST(TestSdeSeirvv, Simulation)
{
    // make a single integration step via a simulation
    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<ScalarType>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(9)) // one call for each flow
        .WillRepeatedly(testing::Return(.0)); // "disable" noise term, as it only adds sqrts of flows

    auto sim = mio::StochasticSimulation(sseirvv_testing_model(), 0.0, 1.0);
    sim.advance(1);

    EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // Stores initial value and single step.

    auto expected_result = Eigen::Vector<double, 10>{0, 1, 1.25, 0.75, 1.25, 0.75, 1.5, 1.25, 0.75, 1.5};
    EXPECT_EQ(sim.get_result().get_last_value(), expected_result);
}

TEST(TestSdeSeirvv, check_constraints_parameters)
{
    // check parameters.check_constraints
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
    // check parameters.apply_constraints
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