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
#include "sde_sirs/model.h"
#include "memilio/compartments/stochastic_simulation.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

const mio::ssirs::Model<double>& ssirs_testing_model()
{
    static mio::ssirs::Model<double> model;
    model.populations.array().setConstant(1);
    model.parameters.set<mio::ssirs::TimeImmune<double>>(4);
    {
        model.parameters.set<mio::ssirs::TimeInfected<double>>(2);
        model.parameters.set<mio::ssirs::TransmissionProbabilityOnContact<double>>(1);
        model.parameters.get<mio::ssirs::ContactPatterns<double>>().get_baseline()(0, 0) = 3;
    }
    return model;
}

TEST(TestSdeSirs, Model)
{
    // check get_flows and get_noise
    const Eigen::Vector3d y = Eigen::Vector3d::Constant(1);
    Eigen::Vector3d flows   = Eigen::Vector3d::Constant(1);
    Eigen::Vector3d noise   = Eigen::Vector3d::Constant(1);

    ssirs_testing_model().get_flows(y, y, 0, flows);
    auto expected_flows = Eigen::Vector3d{1, 0.5, 0.25};
    EXPECT_EQ(flows, expected_flows);

    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<double>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(3)) // one call for each flow
        .WillRepeatedly(testing::Return(0.5));

    ssirs_testing_model().get_noise(y, y, 0, noise);
    Eigen::Vector3d expected_noise;
    ssirs_testing_model().get_derivatives(expected_flows.array().sqrt() * 0.5, expected_noise);
    EXPECT_EQ(noise, expected_noise);
}

TEST(TestSdeSirs, Simulation)
{
    // make a single integration step via a simulation
    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<double>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(3)) // one call for each flow
        .WillRepeatedly(testing::Return(.0)); // "disable" noise term, as it only adds sqrts of flows

    auto sim = mio::StochasticSimulation<double, mio::ssirs::Model<double>>(ssirs_testing_model(), 0.0, 1.0);
    sim.advance(1);

    EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // stores initial value and single step

    auto expected_result = Eigen::Vector3d{0.25, 1.5, 1.25};
    EXPECT_EQ(sim.get_result().get_last_value(), expected_result);
}

TEST(TestSdeSirs, check_constraints_parameters)
{
    // check parameters.check_constraints
    mio::ssirs::Model<double>::ParameterSet parameters;
    parameters.set<mio::ssirs::TimeInfected<double>>(6);
    parameters.set<mio::ssirs::TimeImmune<double>>(6);
    parameters.set<mio::ssirs::TransmissionProbabilityOnContact<double>>(0.04);
    parameters.get<mio::ssirs::ContactPatterns<double>>().get_baseline()(0, 0) = 10;
    parameters.set<mio::ssirs::StartDay<double>>(30);
    parameters.set<mio::ssirs::Seasonality<double>>(0.3);

    // model.check_constraints() combines the functions from population and parameters.
    // We only want to test the functions for the parameters defined in parameters.h
    EXPECT_EQ(parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    parameters.set<mio::ssirs::TimeInfected<double>>(0);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::ssirs::TimeInfected<double>>(6);
    parameters.set<mio::ssirs::TimeImmune<double>>(0);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::ssirs::TimeImmune<double>>(6);
    parameters.set<mio::ssirs::TransmissionProbabilityOnContact<double>>(10.);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::ssirs::TransmissionProbabilityOnContact<double>>(0.04);
    parameters.set<mio::ssirs::Seasonality<double>>(-2.);
    EXPECT_EQ(parameters.check_constraints(), 1);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestSdeSirs, apply_constraints_parameters)
{
    // check parameters.apply_constraints
    const double tol_times = 1e-1;
    mio::ssirs::Model<double>::ParameterSet parameters;
    parameters.set<mio::ssirs::TimeInfected<double>>(6);
    parameters.set<mio::ssirs::TimeImmune<double>>(6);
    parameters.set<mio::ssirs::TransmissionProbabilityOnContact<double>>(0.04);
    parameters.get<mio::ssirs::ContactPatterns<double>>().get_baseline()(0, 0) = 10;
    parameters.set<mio::ssirs::StartDay<double>>(30);
    parameters.set<mio::ssirs::Seasonality<double>>(0.3);

    EXPECT_EQ(parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    parameters.set<mio::ssirs::TimeInfected<double>>(-2.5);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_EQ(parameters.get<mio::ssirs::TimeInfected<double>>(), tol_times);

    parameters.set<mio::ssirs::TimeImmune<double>>(-2.5);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_EQ(parameters.get<mio::ssirs::TimeImmune<double>>(), tol_times);

    parameters.set<mio::ssirs::TransmissionProbabilityOnContact<double>>(10.);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_NEAR(parameters.get<mio::ssirs::TransmissionProbabilityOnContact<double>>(), 0.0, 1e-14);

    parameters.set<mio::ssirs::Seasonality<double>>(-2.);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_NEAR(parameters.get<mio::ssirs::Seasonality<double>>(), 0.0, 1e-14);
    mio::set_log_level(mio::LogLevel::warn);
}
