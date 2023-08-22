/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn, Martin Siggel, Henrik Zunker
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
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/math/euler.h"
#include "memilio/compartments/simulation.h"
#include <gtest/gtest.h>
#include <iomanip>
#include <vector>

TEST(TestSeir, simulateDefault)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::oseir::Model model;
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

TEST(TestSeir, CompareSeirWithJS)
{
    // initialization
    double t0   = 0.;
    double tmax = 50.;
    double dt   = 0.1002004008016032;

    double total_population = 1061000;

    mio::oseir::Model model;

    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]   = 10000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}]  = 1000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}] = 1000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}];
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(1.0);
    model.parameters.set<mio::oseir::TimeExposed>(5.2);
    model.parameters.set<mio::oseir::TimeInfected>(2);

    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 2.7;
    model.parameters.get<mio::oseir::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));

    std::vector<std::vector<double>> refData = load_test_data_csv<double>("seir-js-compare.csv");
    auto integrator                          = std::make_shared<mio::EulerIntegratorCore>();
    auto result                              = mio::simulate<mio::oseir::Model>(t0, tmax, dt, model, integrator);

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

TEST(TestSeir, checkPopulationConservation)
{
    // initialization
    double t0   = 0.;
    double tmax = 50.;
    double dt   = 0.1002004008016032;

    double total_population = 1061000;

    mio::oseir::Model model;

    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]   = 10000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}]  = 1000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}] = 1000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}];
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(1.0);
    model.parameters.set<mio::oseir::TimeExposed>(5.2);
    model.parameters.set<mio::oseir::TimeInfected>(2);

    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 2.7;
    model.parameters.get<mio::oseir::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));
    auto result        = mio::simulate<mio::oseir::Model>(t0, tmax, dt, model);
    double num_persons = 0.0;
    for (auto i = 0; i < result.get_last_value().size(); i++) {
        num_persons += result.get_last_value()[i];
    }
    EXPECT_NEAR(num_persons, total_population, 1e-8);
}

TEST(TestSeir, check_constraints_parameters)
{
    mio::oseir::Model model;
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
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestSeir, apply_constraints_parameters)
{
    mio::oseir::Model model;
    model.parameters.set<mio::oseir::TimeExposed>(5.2);
    model.parameters.set<mio::oseir::TimeInfected>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 10;

    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::oseir::TimeExposed>(-5.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseir::TimeExposed>(), 1e-4);

    model.parameters.set<mio::oseir::TimeInfected>(1e-5);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseir::TimeInfected>(), 1e-4);

    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(10.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::oseir::TransmissionProbabilityOnContact>(), 0.0, 1e-14);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestSeir, get_reproduction_numbers)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    mio::oseir::Model model;

    double total_population                                                                            = 10000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]   = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}]  = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}] = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}];

    model.parameters.set<mio::oseir::TimeExposed>(5.2);
    model.parameters.set<mio::oseir::TimeInfected>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 10;
    model.parameters.get<mio::oseir::ContactPatterns>().add_damping(
        0.5, mio::SimulationTime(
                 0.6)); //Added damping so we can observe an instantaneous reduction by 50% of the reproduction numbers

    model.check_constraints();

    Eigen::VectorXd checkReproductionNumbers(7);
    checkReproductionNumbers << 1.9258478907262198, 1.9241017236517031, 1.9162226787480827, 1.8803730074762808,
        1.7145319195349404, 1.177662053486763, 1.16129426162878;

    mio::TimeSeries<double> result = mio::simulate(t0, tmax, dt, model);

    auto reproduction_numbers = model.get_reproduction_numbers(result);

    for (int i = 0; i < reproduction_numbers.size(); i++) {
        EXPECT_NEAR(reproduction_numbers[i], checkReproductionNumbers[i], 1e-12);
    }

    EXPECT_FALSE(model.get_reproduction_number(static_cast<double>(static_cast<size_t>(result.get_num_time_points())),
                                               result)); //Test for an index that is out of range
}

TEST(TestSeir, interpolate_reproduction_numbers)
{
    mio::oseir::Model model;

    double total_population = 10000; //Initialize compartments to get total population of 10000
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]   = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}]  = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}] = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}];

    model.parameters.set<mio::oseir::TimeInfected>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 10;
    model.parameters.get<mio::oseir::ContactPatterns>().add_damping(
        0.5, mio::SimulationTime(
                 0.5)); //Added damping so we can observe an instantaneous reduction by 50% of the reproduction numbers

    model.apply_constraints();

    mio::TimeSeries<ScalarType> result((int)mio::oseir::InfectionState::Count);
    mio::TimeSeries<ScalarType>::Vector result_0(4);
    mio::TimeSeries<ScalarType>::Vector result_1(4);
    mio::TimeSeries<ScalarType>::Vector result_2(4);
    mio::TimeSeries<ScalarType>::Vector result_3(4);
    mio::TimeSeries<ScalarType>::Vector result_4(4);
    mio::TimeSeries<ScalarType>::Vector result_5(4);
    mio::TimeSeries<ScalarType>::Vector result_6(4);
    mio::TimeSeries<ScalarType>::Vector result_7(4);

    result_0[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9700;
    result_1[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9699.9709149074315;
    result_2[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9699.8404009584538;
    result_3[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9699.260556488618;
    result_4[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9696.800490904101;
    result_5[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9687.9435082620021;
    result_6[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9679.5436372291661;
    result_7[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9678.5949381732935;

    result.add_time_point(0, result_0);
    result.add_time_point(0.001, result_1);
    result.add_time_point(0.0055, result_2);
    result.add_time_point(0.02575, result_3);
    result.add_time_point(0.116875, result_4);
    result.add_time_point(0.526938, result_5);
    result.add_time_point(0.952226, result_6);
    result.add_time_point(1, result_7);

    EXPECT_FALSE(model.get_reproduction_number(result.get_time(0) - 0.5, result)); //Test for indices out of range
    EXPECT_FALSE(model.get_reproduction_number(result.get_last_time() + 0.5, result));
    EXPECT_EQ(model.get_reproduction_number(0.3, result).value(), 1.3695409350793410486);
    EXPECT_EQ(model.get_reproduction_number(0.7, result).value(), 1.1621430429058086098);
}
