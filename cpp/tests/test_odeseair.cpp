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
#include "load_test_data.h"
#include <gtest/gtest.h>

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

    const double N = 327167434; // total population of the United States

    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] =
        0.9977558755803503;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}] =
        0.0003451395725394549;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] =
        0.00037846880968213874;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}] =
        (337072.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] =
        (17448.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Perished)}] = (9619.0 / N);
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
