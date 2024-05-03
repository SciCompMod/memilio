/*
* Copyright (C) 2020-2024 MEmilio
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

        total_population = 51; // total population of the United States

        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}]  = 20;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}]      = 1;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] = 2;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}]     = 4;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}]    = 8;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Perished)}]     = 16;
        model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::ObjectiveFunction)}] =
            0.0;

        model.parameters.set<mio::oseair::AlphaA<double>>(0.5);
        model.parameters.set<mio::oseair::AlphaI<double>>(0.5);
        model.parameters.get<mio::oseair::Kappa<double>>()          = 0.5;
        model.parameters.get<mio::oseair::Beta<double>>()           = 0.5;
        model.parameters.get<mio::oseair::Mu<double>>()             = 0.5;
        model.parameters.get<mio::oseair::TLatentInverse<double>>() = 0.5;
        model.parameters.get<mio::oseair::Rho<double>>()            = 0.5;
        model.parameters.get<mio::oseair::Gamma<double>>()          = 0.5;
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
    auto result = mio::simulate<double, mio::oseair::Model<double>>(t0, 50, dt, model);
    double num_persons =
        result.get_last_value().sum() - result.get_last_value()[(size_t)mio::oseair::InfectionState::ObjectiveFunction];
    EXPECT_NEAR(num_persons, total_population, 1e-8);
}

TEST_F(ModelTestOdeSeair, get_derivatives)
{
    auto dydt_default = Eigen::VectorXd(Eigen::Index(mio::oseair::InfectionState::Count));
    dydt_default.setZero();
    auto y0 = model.get_initial_values();
    model.get_derivatives(y0, y0, 0, dydt_default);

    EXPECT_NEAR(dydt_default[0], -56, 1e-12);
    EXPECT_NEAR(dydt_default[1], 59.5, 1e-12);
    EXPECT_NEAR(dydt_default[2], -3, 1e-12);
    EXPECT_NEAR(dydt_default[3], -1, 1e-12);
    EXPECT_NEAR(dydt_default[4], -1.5, 1e-12);
    EXPECT_NEAR(dydt_default[5], 2, 1e-12);
    EXPECT_NEAR(dydt_default[6], 0.05, 1e-12);
}

TEST_F(ModelTestOdeSeair, Simulation)
{
    dt                                                      = 1;
    std::shared_ptr<mio::IntegratorCore<double>> integrator = std::make_shared<mio::EulerIntegratorCore<double>>();
    auto sim                                                = simulate(t0, tmax, dt, model, integrator);

    EXPECT_EQ(sim.get_num_time_points(), 2);

    const auto& results_t1 = sim.get_last_value();
    auto results_t1_vec    = std::vector<double>(results_t1.data(), results_t1.data() + results_t1.size());
    EXPECT_NEAR(results_t1[0], -36, 1e-12);
    EXPECT_NEAR(results_t1[1], 60.5, 1e-12);
    EXPECT_NEAR(results_t1[2], 1, 1e-12);
    EXPECT_NEAR(results_t1[3], 7, 1e-12);
    EXPECT_NEAR(results_t1[4], 0.5, 1e-12);
    EXPECT_NEAR(results_t1[5], 18, 1e-12);
    EXPECT_NEAR(results_t1[6], 0.05, 1e-12);
}

TEST(TestOdeSeair, compareWithPreviousRun)
{
    double t0   = 0;
    double tmax = 10;
    double dt   = 0.3;

    mio::oseair::Model<double> model;

    // set parameters to default values just for code covereage

    model.parameters.set<mio::oseair::AlphaA<double>>(0.2);
    model.parameters.set<mio::oseair::AlphaI<double>>(0.2);

    model.parameters.get<mio::oseair::Kappa<double>>()          = 0.2;
    model.parameters.get<mio::oseair::Beta<double>>()           = 0.0067;
    model.parameters.get<mio::oseair::Mu<double>>()             = 0.0041;
    model.parameters.get<mio::oseair::TLatentInverse<double>>() = 0.5;
    model.parameters.get<mio::oseair::Rho<double>>()            = 0.1;
    model.parameters.get<mio::oseair::Gamma<double>>()          = 0.0;

    // set initial values

    const double total_population = 327167434; // total population of the United States

    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] =
        0.9977558755803503;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}] =
        0.0003451395725394549;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] =
        0.00037846880968213874;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}] =
        (337072.0 / total_population);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] =
        (17448.0 / total_population);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Perished)}] =
        (9619.0 / total_population);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::ObjectiveFunction)}] = 0.0;

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
