/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Ralf Hannemann-Tamas
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

#include "ode_seair/model.h"
#include "ode_seair/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/euler.h"
#include "load_test_data.h"
#include <gtest/gtest.h>

class ModelTestOdeSeair : public testing::Test
{

public:
    double t0;
    double tmax;
    double dt;
    double total_population;
    mio::oseair::Model<double> model;

protected:
    void SetUp() override
    {
        t0   = 0.;
        tmax = 1.;
        dt   = 0.1;

        total_population                                                                                        = 100;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}]      = 10;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] = 20;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}]     = 5;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}]    = 10;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Dead)}]         = 10;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] =
            total_population -
            model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}] -
            model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] -
            model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}] -
            model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] -
            model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Dead)}];

        model.parameters.set<mio::oseair::SocialDistancing<double>>(.5);
        model.parameters.set<mio::oseair::Quarantined<double>>(.5);
        model.parameters.get<mio::oseair::TestingRate<double>>()                  = 1.;
        model.parameters.get<mio::oseair::RecoveryRate<double>>()                 = 1.;
        model.parameters.get<mio::oseair::DeathRate<double>>()                    = 1.;
        model.parameters.get<mio::oseair::TimeExposed<double>>()                  = 1.;
        model.parameters.get<mio::oseair::RecoveryRateFromAsymptomatic<double>>() = 1.;
        model.parameters.get<mio::oseair::TimeRecoveredInv<double>>()             = 1.;
    }
};

TEST_F(ModelTestOdeSeair, simulateDefault)
{
    EXPECT_EQ(model.parameters.check_constraints(), 0);
    auto result = mio::simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST_F(ModelTestOdeSeair, checkPopulationConservation)
{
    auto result        = mio::simulate<double, mio::oseair::Model<double>>(t0, 50, dt, model);
    double num_persons = result.get_last_value().sum();
    EXPECT_NEAR(num_persons, total_population, 1e-8);
}

TEST_F(ModelTestOdeSeair, get_derivatives)
{
    auto dydt_default = Eigen::VectorXd(Eigen::Index(mio::oseair::InfectionState::Count));
    dydt_default.setZero();
    auto y0 = model.get_initial_values();
    model.get_derivatives(y0, y0, 0, dydt_default);

    EXPECT_NEAR(dydt_default[(size_t)mio::oseair::InfectionState::Susceptible], 4.375, 1e-12);
    EXPECT_NEAR(dydt_default[(size_t)mio::oseair::InfectionState::Exposed], -4.375, 1e-12);
    EXPECT_NEAR(dydt_default[(size_t)mio::oseair::InfectionState::Asymptomatic], -30, 1e-12);
    EXPECT_NEAR(dydt_default[(size_t)mio::oseair::InfectionState::Infected], 10, 1e-12);
    EXPECT_NEAR(dydt_default[(size_t)mio::oseair::InfectionState::Recovered], 15, 1e-12);
    EXPECT_NEAR(dydt_default[(size_t)mio::oseair::InfectionState::Dead], 5, 1e-12);
}

TEST_F(ModelTestOdeSeair, Simulation)
{
    dt                                                      = 1;
    std::shared_ptr<mio::IntegratorCore<double>> integrator = std::make_shared<mio::EulerIntegratorCore<double>>();
    auto sim                                                = simulate(t0, tmax, dt, model, integrator);

    EXPECT_EQ(sim.get_num_time_points(), 2);

    const auto& results_t1 = sim.get_last_value();
    EXPECT_NEAR(results_t1[(size_t)mio::oseair::InfectionState::Susceptible], 49.375, 1e-12);
    EXPECT_NEAR(results_t1[(size_t)mio::oseair::InfectionState::Exposed], 5.625, 1e-12);
    EXPECT_NEAR(results_t1[(size_t)mio::oseair::InfectionState::Asymptomatic], -10, 1e-12);
    EXPECT_NEAR(results_t1[(size_t)mio::oseair::InfectionState::Infected], 15, 1e-12);
    EXPECT_NEAR(results_t1[(size_t)mio::oseair::InfectionState::Recovered], 25, 1e-12);
    EXPECT_NEAR(results_t1[(size_t)mio::oseair::InfectionState::Dead], 15, 1e-12);
}

TEST(TestOdeSeair, compareWithPreviousRun)
{
    double t0   = 0;
    double tmax = 10;
    double dt   = 0.3;

    mio::oseair::Model<double> model;

    // Set the parameters to the values used.
    model.parameters.set<mio::oseair::SocialDistancing<double>>(0.2);
    model.parameters.set<mio::oseair::Quarantined<double>>(0.2);

    model.parameters.get<mio::oseair::TestingRate<double>>()                  = 0.2;
    model.parameters.get<mio::oseair::RecoveryRate<double>>()                 = 0.0067;
    model.parameters.get<mio::oseair::DeathRate<double>>()                    = 0.0041;
    model.parameters.get<mio::oseair::TimeExposed<double>>()                  = 2.;
    model.parameters.get<mio::oseair::RecoveryRateFromAsymptomatic<double>>() = 0.1;
    model.parameters.get<mio::oseair::TimeRecoveredInv<double>>()             = 0.0;

    // set initial values
    const double total_population = 327167434; // total population of the United States
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}] =
        0.0003451395725394549;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] =
        0.00037846880968213874;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}] =
        (337072.0 / total_population);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] =
        (17448.0 / total_population);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Dead)}] =
        (9619.0 / total_population);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] =
        1 - model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}] -
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] -
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] -
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Dead)}];

    // check constraints
    model.apply_constraints();

    // simulate
    mio::TimeSeries<double> seair = mio::simulate<double, mio::oseair::Model<double>>(t0, tmax, dt, model);

    auto compare = load_test_data_csv<double>("seair-compare.csv");

    ASSERT_EQ(compare.size(), static_cast<size_t>(seair.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(seair.get_num_elements()) + 1) << "at row " << i;
        EXPECT_NEAR(seair.get_time(i), compare[i][0], 1e-10) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            EXPECT_NEAR(seair.get_value(i)[j - 1], compare[i][j], 1e-10) << " at row " << i;
        }
    }
}

// Test if check_constraints() work as expected.
TEST(TestOdeSeair, testParameterConstraints)
{ // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    mio::oseair::Parameters<> parameters;
    double negative_value = -0.5;

    // Check if default parameters are valid.
    bool constraint_check = parameters.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Check improper SocialDistancing.
    parameters.get<mio::oseair::SocialDistancing<>>() = negative_value;
    constraint_check                                  = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters.get<mio::oseair::SocialDistancing<>>() = 0.5;

    // Check improper Quarantined.
    parameters.get<mio::oseair::Quarantined<>>() = negative_value;
    constraint_check                             = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters.get<mio::oseair::Quarantined<>>() = 0.5;

    // Check improper TimeExposed.
    parameters.get<mio::oseair::TimeExposed<>>() = 0;
    constraint_check                             = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters.get<mio::oseair::TimeExposed<>>() = 3.1;

    // Check improper RecoveryRateFromAsymptomatic.
    parameters.get<mio::oseair::RecoveryRateFromAsymptomatic<>>() = negative_value;
    constraint_check                                              = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters.get<mio::oseair::RecoveryRateFromAsymptomatic<>>() = 0.7 / 5.;

    // Check improper TestingRate.
    parameters.get<mio::oseair::TestingRate<>>() = negative_value;
    constraint_check                             = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters.get<mio::oseair::TestingRate<>>() = 0.3 / 5.;

    // Check improper RecoveryRate.
    parameters.get<mio::oseair::RecoveryRate<>>() = negative_value;
    constraint_check                              = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters.get<mio::oseair::RecoveryRate<>>() = 0.5 / 7.;

    // Check improper DeathRate.
    parameters.get<mio::oseair::DeathRate<>>() = negative_value;
    constraint_check                           = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters.get<mio::oseair::DeathRate<>>() = 0.5 / 7.;

    // Check improper TimeRecoveredInv.
    parameters.get<mio::oseair::TimeRecoveredInv<>>() = negative_value;
    constraint_check                                  = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters.get<mio::oseair::TimeRecoveredInv<>>() = 0.;

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);

    // Check with correct parameters.
    constraint_check = parameters.check_constraints();
    EXPECT_FALSE(constraint_check);
}
