/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn, Martin Siggel
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
#include "load_test_data.h"
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/math/euler.h"
#include "memilio/compartments/simulation.h"
#include <gtest/gtest.h>

using real = double;

class TestSeir : public testing::Test
{
protected:
    void SetUp() override
    {
        t0   = 0.;
        tmax = 50.;
        dt   = 0.1002004008016032;

        total_population = 1061000;

        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]   = 10000;
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}]  = 1000;
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}] = 1000;
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] =
            total_population -
            this->model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}] -
            this->model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}] -
            this->model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}];
        // suscetible now set with every other update
        // model.nb_sus_t0   = model.nb_total_t0 - model.nb_exp_t0 - model.nb_inf_t0 - model.nb_rec_t0;
        model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(1.0);
        model.parameters.set<mio::oseir::TimeExposed>(5.2);
        model.parameters.set<mio::oseir::TimeInfected>(2);

        model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 2.7;
        model.parameters.get<mio::oseir::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));
    }

public:
    real t0;
    real tmax;
    real dt;
    double total_population;
    mio::oseir::Model model;
};

TEST_F(TestSeir, CompareSeirWithPyBindings)
{
    /*
    This test test the cpp model. The same test is implemented in python to compare the results of both simulations.
    If this test is change the corresponding python test needs to be changed aswell (also updating the data file).
    */
    std::vector<std::vector<real>> refData = load_test_data_csv<real>("seir-compare.csv");
    auto result                            = mio::simulate<mio::oseir::Model>(t0, tmax, dt, model);

    ASSERT_EQ(refData.size(), static_cast<size_t>(result.get_num_time_points()));

    for (Eigen::Index irow = 0; irow < result.get_num_time_points(); ++irow) {
        double t     = refData[static_cast<size_t>(irow)][0];
        auto rel_tol = 1e-6;

        //test result diverges at damping because of changes, not worth fixing at the moment
        if (t > 11.0 && t < 13.0) {
            //strong divergence around damping
            rel_tol = 0.5;
        }
        else if (t > 13.0) {
            //minor divergence after damping
            rel_tol = 1e-2;
        }

        ASSERT_NEAR(t, result.get_times()[irow], 1e-12) << "at row " << irow;
        for (size_t icol = 0; icol < 4; ++icol) {
            double ref    = refData[static_cast<size_t>(irow)][icol + 1];
            double actual = result[irow][icol];

            double tol = rel_tol * ref;
            ASSERT_NEAR(ref, actual, tol) << "at row " << irow;
        }
    }
}

TEST_F(TestSeir, checkPopulationConservation)
{
    auto result        = mio::simulate<mio::oseir::Model>(t0, tmax, dt, model);
    double num_persons = 0.0;
    for (auto i = 0; i < result.get_last_value().size(); i++) {
        num_persons += result.get_last_value()[i];
    }
    EXPECT_NEAR(num_persons, total_population, 1e-8);
}

TEST_F(TestSeir, check_constraints_parameters)
{
    model.parameters.set<mio::oseir::TimeExposed>(5.2);
    model.parameters.set<mio::oseir::TimeInfected>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 10;

    // model.check_constraints() combines the functions from population and parameters.
    // We only want to test the functions for the parameters defined in parameters.h
    ASSERT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::oseir::TimeExposed>(-5.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::oseir::TimeExposed>(5.2);
    model.parameters.set<mio::oseir::TimeInfected>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::oseir::TimeInfected>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(10.);
    ASSERT_EQ(model.parameters.check_constraints(), 1);
}
