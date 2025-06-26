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
#include "memilio/compartments/stochastic_simulation.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

const mio::ssir::Model& ssir_testing_model()
{
    static mio::ssir::Model model;
    model.populations.array().setConstant(1);
    {
        model.parameters.set<mio::ssir::TimeInfected>(0.5);
        model.parameters.set<mio::ssir::TransmissionProbabilityOnContact>(1);
        model.parameters.get<mio::ssir::ContactPatterns>().get_baseline()(0, 0) = 3;
    }
    return model;
}

TEST(TestSdeSir, Model)
{
    // check get_flows and get_noise
    const Eigen::Vector3d y = Eigen::Vector3d::Constant(1);
    Eigen::Vector2d rhs     = Eigen::Vector2d::Constant(1);

    ssir_testing_model().get_flows(y, y, 0, rhs);
    auto expected_flows = Eigen::Vector2d{1, 2};
    EXPECT_EQ(rhs, expected_flows);

    ssir_testing_model().get_noise(y, y, 0, rhs);
    auto expected_noise = expected_flows.array().sqrt().matrix().eval();
    EXPECT_EQ(rhs, expected_noise);
}

TEST(TestSdeSir, Simulation)
{
    // make a single integration step via a simulation
    ScopedMockDistribution<
        testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<double>>>>>
        normal_dist_mock;

    EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
        .Times(testing::Exactly(3)) // one call for each compartment
        .WillRepeatedly(testing::Return(.0)); // "disable" noise term, as it only adds sqrts of flows

    auto sim = mio::StochasticSimulation(ssir_testing_model(), 0.0, 1.0);
    sim.advance(1.0);

    EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // stores initial value and single step

    auto expected_result = Eigen::Vector3d{0, 0, 3};
    EXPECT_EQ(sim.get_result().get_last_value(), expected_result);
}

// TEST(TestSdeSir, FlowSimulation)
// {
//     // make a single integration step via a flow simulation
//     ScopedMockDistribution<
//         testing::StrictMock<MockDistribution<mio::DistributionAdapter<std::normal_distribution<double>>>>>
//         normal_dist_mock;

//     EXPECT_CALL(normal_dist_mock.get_mock(), invoke)
//         .Times(testing::Exactly(2)) // one call for each flow
//         .WillRepeatedly(testing::Return(.0)); // "disable" noise term, as it only adds sqrts of flows

//     auto sim = mio::StochasticFlowSimulation(ssir_testing_model(), 0.0, 1.0);
//     sim.advance(1.0);

//     EXPECT_EQ(sim.get_result().get_num_time_points(), 2); // stores initial value and single step

//     auto expected_result = Eigen::Vector3d{0, 0, 3};
//     EXPECT_EQ(sim.get_result().get_last_value(), expected_result);

//     auto expected_flows = Eigen::Vector2d{1, 2};
//     EXPECT_EQ(sim.get_flows().get_last_value(), expected_flows);
// }

TEST(TestSdeSir, check_constraints_parameters)
{
    mio::ssir::Model::ParameterSet parameters;
    parameters.set<mio::ssir::TimeInfected>(6);
    parameters.set<mio::ssir::TransmissionProbabilityOnContact>(0.04);
    parameters.get<mio::ssir::ContactPatterns>().get_baseline()(0, 0) = 10;

    // model.check_constraints() combines the functions from population and parameters.
    // We only want to test the functions for the parameters defined in parameters.h
    EXPECT_EQ(parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    parameters.set<mio::ssir::TimeInfected>(0);
    EXPECT_EQ(parameters.check_constraints(), 1);

    parameters.set<mio::ssir::TimeInfected>(6);
    parameters.set<mio::ssir::TransmissionProbabilityOnContact>(10.);
    EXPECT_EQ(parameters.check_constraints(), 1);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestSdeSir, apply_constraints_parameters)
{
    const double tol_times = 1e-1;
    mio::ssir::Model::ParameterSet parameters;
    parameters.set<mio::ssir::TimeInfected>(6);
    parameters.set<mio::ssir::TransmissionProbabilityOnContact>(0.04);
    parameters.get<mio::ssir::ContactPatterns>().get_baseline()(0, 0) = 10;

    EXPECT_EQ(parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    parameters.set<mio::ssir::TimeInfected>(-2.5);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_EQ(parameters.get<mio::ssir::TimeInfected>(), tol_times);

    parameters.set<mio::ssir::TransmissionProbabilityOnContact>(10.);
    EXPECT_EQ(parameters.apply_constraints(), 1);
    EXPECT_NEAR(parameters.get<mio::ssir::TransmissionProbabilityOnContact>(), 0.0, 1e-14);
    mio::set_log_level(mio::LogLevel::warn);
}
