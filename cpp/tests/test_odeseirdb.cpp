/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
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

#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "ode_seirdb/infection_state.h"
#include "ode_seirdb/model.h"
#include "ode_seirdb/parameters.h"

#include <gtest/gtest.h>

TEST(TestOdeSeirdb, simulateDefault)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::oseirdb::Model<double> model(1);
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

class ModelTestOdeSeirdb : public testing::Test
{
public:
    ModelTestOdeSeirdb()
        : model(1)
    {
    }
    double t0;
    double tmax;
    double dt;
    double total_population;
    mio::oseirdb::Model<double> model;

protected:
    void SetUp() override
    {
        t0   = 0.;
        tmax = 50.;
        dt   = 0.25;

        total_population = 100000.;

        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Exposed}]   = 1000;
        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Infected}]  = 1000;
        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Recovered}] = 1000;
        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Dead}]      = 100;
        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Buried}]    = 80;
        model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Susceptible}] =
            total_population - model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Exposed}] -
            model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Infected}] -
            model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Recovered}] -
            model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Dead}] -
            model.populations[{mio::AgeGroup(0), mio::oseirdb::InfectionState::Buried}];

        model.parameters.set<mio::oseirdb::TransmissionProbabilityOnContact<double>>(0.1);
        model.parameters.set<mio::oseirdb::TransmissionProbabilityFromDead<double>>(0.01);
        model.parameters.set<mio::oseirdb::ProbabilityToRecover<double>>(0.6);
        model.parameters.set<mio::oseirdb::TimeExposed<double>>(5.2);
        model.parameters.set<mio::oseirdb::TimeInfected<double>>(6.0);
        model.parameters.set<mio::oseirdb::TimeToBurial<double>>(4.0);

        mio::ContactMatrixGroup<double>& contact_matrix =
            model.parameters.get<mio::oseirdb::ContactPatterns<double>>().get_cont_freq_mat();
        contact_matrix[0].get_baseline().setConstant(2.0);
        contact_matrix[0].add_damping(0.5, mio::SimulationTime<double>(10.));
    }
};

TEST_F(ModelTestOdeSeirdb, checkPopulationConservation)
{
    auto result        = mio::simulate<double, mio::oseirdb::Model<double>>(t0, tmax, dt, model);
    double num_persons = result.get_last_value().sum();
    EXPECT_NEAR(num_persons, total_population, 1e-6);
}

TEST_F(ModelTestOdeSeirdb, check_constraints_parameters)
{
    ASSERT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);

    model.parameters.set<mio::oseirdb::TimeExposed<double>>(-1.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);
    model.parameters.set<mio::oseirdb::TimeExposed<double>>(5.2);

    model.parameters.set<mio::oseirdb::TimeInfected<double>>(0.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);
    model.parameters.set<mio::oseirdb::TimeInfected<double>>(3.0);

    model.parameters.set<mio::oseirdb::TimeToBurial<double>>(0.0);
    EXPECT_EQ(model.parameters.check_constraints(), 1);
    model.parameters.set<mio::oseirdb::TimeToBurial<double>>(2.0);

    model.parameters.set<mio::oseirdb::TransmissionProbabilityOnContact<double>>(1.5);
    EXPECT_EQ(model.parameters.check_constraints(), 1);
    model.parameters.set<mio::oseirdb::TransmissionProbabilityOnContact<double>>(0.05);

    model.parameters.set<mio::oseirdb::TransmissionProbabilityFromDead<double>>(-0.1);
    EXPECT_EQ(model.parameters.check_constraints(), 1);
    model.parameters.set<mio::oseirdb::TransmissionProbabilityFromDead<double>>(0.02);

    model.parameters.set<mio::oseirdb::ProbabilityToRecover<double>>(1.5);
    EXPECT_EQ(model.parameters.check_constraints(), 1);
    model.parameters.set<mio::oseirdb::ProbabilityToRecover<double>>(0.6);

    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSeirdb, apply_constraints_parameters)
{
    constexpr double tol_times = 1e-1;
    mio::oseirdb::Model<double> model(1);

    model.parameters.set<mio::oseirdb::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseirdb::TimeInfected<double>>(6.0);
    model.parameters.set<mio::oseirdb::TimeToBurial<double>>(2.0);
    model.parameters.set<mio::oseirdb::TransmissionProbabilityOnContact<double>>(0.04);
    model.parameters.set<mio::oseirdb::TransmissionProbabilityFromDead<double>>(0.03);
    model.parameters.set<mio::oseirdb::ProbabilityToRecover<double>>(0.4);

    mio::ContactMatrixGroup<double>& contact_matrix =
        model.parameters.get<mio::oseirdb::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(5.0);

    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::oseirdb::TimeExposed<double>>(-5.0);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirdb::TimeExposed<double>>()[(mio::AgeGroup)0], tol_times);

    model.parameters.set<mio::oseirdb::TimeInfected<double>>(1e-6);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirdb::TimeInfected<double>>()[(mio::AgeGroup)0], tol_times);

    model.parameters.set<mio::oseirdb::TimeToBurial<double>>(1e-7);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseirdb::TimeToBurial<double>>()[(mio::AgeGroup)0], tol_times);

    model.parameters.set<mio::oseirdb::TransmissionProbabilityOnContact<double>>(10.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_DOUBLE_EQ(model.parameters.get<mio::oseirdb::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0],
                     0.0);

    model.parameters.set<mio::oseirdb::TransmissionProbabilityFromDead<double>>(-1.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_DOUBLE_EQ(model.parameters.get<mio::oseirdb::TransmissionProbabilityFromDead<double>>()[(mio::AgeGroup)0],
                     0.0);

    model.parameters.set<mio::oseirdb::ProbabilityToRecover<double>>(-0.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_DOUBLE_EQ(model.parameters.get<mio::oseirdb::ProbabilityToRecover<double>>()[(mio::AgeGroup)0], 0.0);
    mio::set_log_level(mio::LogLevel::warn);
}
