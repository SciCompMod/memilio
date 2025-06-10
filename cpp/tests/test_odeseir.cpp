/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/math/integrator.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/math/euler.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"

#include <gtest/gtest.h>

#include <memory>

TEST(TestOdeSeir, simulateDefault)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.1;

    mio::oseir::Model<double> model(1);
    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
}

class ModelTestOdeSeir : public testing::Test
{

public:
    ModelTestOdeSeir()
        : model(1)
    {
    }
    double t0;
    double tmax;
    double dt;
    double total_population;
    mio::oseir::Model<double> model;

protected:
    void SetUp() override
    {
        t0   = 0.;
        tmax = 50.;
        dt   = 0.1002004008016032;

        total_population = 1061000;

        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 10000;
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 1000;
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 1000;
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
            total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
            model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
            model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];
        model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(1.0);
        model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
        model.parameters.set<mio::oseir::TimeInfected<double>>(2);

        mio::ContactMatrixGroup& contact_matrix =
            model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
        contact_matrix[0].get_baseline().setConstant(2.7);
        contact_matrix[0].add_damping(0.6, mio::SimulationTime(12.5));
    }
};

TEST_F(ModelTestOdeSeir, checkPopulationConservation)
{

    auto result = mio::simulate<double, mio::oseir::Model<double>>(t0, tmax, dt, model);

    double num_persons = result.get_last_value().sum();
    EXPECT_NEAR(num_persons, total_population, 1e-8);
}

TEST_F(ModelTestOdeSeir, check_constraints_parameters)
{
    model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.04);

    model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat()[0].get_baseline()(0, 0) = 10;

    // model.check_constraints() combines the functions from population and parameters.
    // We only want to test the functions for the parameters defined in parameters.h
    ASSERT_EQ(model.parameters.check_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::oseir::TimeExposed<double>>(-5.2);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(0);
    ASSERT_EQ(model.parameters.check_constraints(), 1);

    model.parameters.set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(10.);
    ASSERT_EQ(model.parameters.check_constraints(), 1);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSeir, apply_constraints_parameters)
{
    const double tol_times = 1e-1;
    mio::oseir::Model<double> model(1);
    model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.04);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(10);

    EXPECT_EQ(model.parameters.apply_constraints(), 0);

    mio::set_log_level(mio::LogLevel::off);
    model.parameters.set<mio::oseir::TimeExposed<double>>(-5.2);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseir::TimeExposed<double>>()[(mio::AgeGroup)0], tol_times);

    model.parameters.set<mio::oseir::TimeInfected<double>>(1e-5);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_EQ(model.parameters.get<mio::oseir::TimeInfected<double>>()[(mio::AgeGroup)0], tol_times);

    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(10.);
    EXPECT_EQ(model.parameters.apply_constraints(), 1);
    EXPECT_NEAR(model.parameters.get<mio::oseir::TransmissionProbabilityOnContact<double>>()[(mio::AgeGroup)0], 0.0,
                1e-14);
    mio::set_log_level(mio::LogLevel::warn);
}

TEST(TestOdeSeir, get_reproduction_numbers)
{
    mio::oseir::Model<double> model(1);

    double total_population                                                      = 10000;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];

    model.parameters.set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.04);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(10);

    model.apply_constraints();

    Eigen::VectorXd checkReproductionNumbers(7);
    checkReproductionNumbers << 2.3280000000000002913, 2.3279906878991880603, 2.3279487809434575851,
        2.3277601483151548756, 2.3269102025388899158, 2.3230580052413736247, 2.3185400624683065729;

    Eigen::VectorXd checkReproductionNumbers2(7);
    checkReproductionNumbers2 << 2.0952000000000001734, 2.0951916191092689878, 2.0951539028491117378,
        2.0949841334836394324, 2.0942191822850007021, 2.0907522047172362178, 2.086686056221475738;

    Eigen::VectorXd checkReproductionNumbers3(7);
    checkReproductionNumbers3 << 1.8623999999999998334, 1.8623925503193501374, 1.8623590247547658905,
        1.8622081186521235452, 1.8615281620311117106, 1.8584464041930985889, 1.854832049974644903;

    mio::TimeSeries<ScalarType> result((int)mio::oseir::InfectionState::Count);
    mio::TimeSeries<ScalarType>::Vector result_0(4);
    mio::TimeSeries<ScalarType>::Vector result_1(4);
    mio::TimeSeries<ScalarType>::Vector result_2(4);
    mio::TimeSeries<ScalarType>::Vector result_3(4);
    mio::TimeSeries<ScalarType>::Vector result_4(4);
    mio::TimeSeries<ScalarType>::Vector result_5(4);
    mio::TimeSeries<ScalarType>::Vector result_6(4);

    result_0.setZero();
    result_1.setZero();
    result_2.setZero();
    result_3.setZero();
    result_4.setZero();
    result_5.setZero();
    result_6.setZero();

    result_0[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9700;
    result_1[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9699.9611995799496071;
    result_2[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9699.7865872644051706;
    result_3[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9699.0006179798110679;
    result_4[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9695.4591772453732119;
    result_5[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9679.4083551723888377;
    result_6[(Eigen::Index)mio::oseir::InfectionState::Susceptible] = 9660.5835936179428245;

    result.add_time_point(0, result_0);
    result.add_time_point(0.0010000000000000000208, result_1);
    result.add_time_point(0.0055000000000000005482, result_2);
    result.add_time_point(0.025750000000000005523, result_3);
    result.add_time_point(0.11687500000000002054, result_4);
    result.add_time_point(0.52693750000000005862, result_5);
    result.add_time_point(1, result_6);

    auto reproduction_numbers = model.get_reproduction_numbers(result);

    for (int i = 0; i < reproduction_numbers.size(); i++) {
        EXPECT_NEAR(reproduction_numbers[i], checkReproductionNumbers[i], 1e-12);
    }

    contact_matrix[0].get_baseline().setConstant(9);

    auto reproduction_numbers2 = model.get_reproduction_numbers(result);

    for (int i = 0; i < reproduction_numbers2.size(); i++) {
        EXPECT_NEAR(reproduction_numbers2[i], checkReproductionNumbers2[i], 1e-12);
    }

    contact_matrix[0].get_baseline().setConstant(8);

    auto reproduction_numbers3 = model.get_reproduction_numbers(result);

    for (int i = 0; i < reproduction_numbers2.size(); i++) {
        EXPECT_NEAR(reproduction_numbers3[i], checkReproductionNumbers3[i], 1e-12);
    }

    EXPECT_FALSE(model.get_reproduction_number(static_cast<double>(static_cast<size_t>(result.get_num_time_points())),
                                               result)); //Test for an index that is out of range
}

TEST(TestOdeSeir, get_reproduction_number)
{
    mio::oseir::Model<double> model(1);

    double total_population = 10000; //Initialize compartments to get total population of 10000
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];

    model.parameters.set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.04);

    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(10);

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
    EXPECT_FALSE(model.get_reproduction_number((size_t)result.get_num_time_points(), result));

    EXPECT_EQ(model.get_reproduction_number((size_t)0, result).value(),
              model.get_reproduction_number(0.0, result).value());

    EXPECT_NEAR(model.get_reproduction_number(0.3, result).value(), 2.3262828383474389859, 1e-12);
    EXPECT_NEAR(model.get_reproduction_number(0.7, result).value(), 2.3242860858116172196, 1e-12);
    EXPECT_NEAR(model.get_reproduction_number(0.0, result).value(), 2.3280000000000002913, 1e-12);

    contact_matrix[0].get_baseline().setConstant(9);
    EXPECT_NEAR(model.get_reproduction_number(0.1, result).value(), 2.0946073086586665113, 1e-12);
    EXPECT_NEAR(model.get_reproduction_number(0.3, result).value(), 2.0936545545126947765, 1e-12);

    contact_matrix[0].get_baseline().setConstant(8);
    EXPECT_NEAR(model.get_reproduction_number(0.2, result).value(), 1.8614409729718137676, 1e-12);
    EXPECT_NEAR(model.get_reproduction_number(0.9, result).value(), 1.858670429549998504, 1e-12);
}

// Test model initialization with total population of 0 and ensure get_flows returns no NaN values
TEST(TestSeir, population_zero_no_nan)
{
    // initialize simple model with total population 0
    mio::oseir::Model<double> model(1);
    model.populations.set_total(0.0);

    // call the get_flows function
    auto dydt_default = Eigen::VectorXd(3);
    dydt_default.setZero();
    auto y0 = model.get_initial_values();
    model.get_flows(y0, y0, 0, dydt_default);

    // check that there are now NaN values in dydt_default
    for (int i = 0; i < dydt_default.size(); i++) {
        EXPECT_FALSE(std::isnan(dydt_default[i]));
    }
}

TEST(TestSeir, get_flows)
{
    mio::oseir::Model<double> model(1);

    constexpr double total_population = 400;

    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];

    model.parameters.set<mio::oseir::TimeExposed<double>>(2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(4);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(1);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(1);
    model.check_constraints();

    auto dydt_default = Eigen::VectorXd(3);
    dydt_default.setZero();
    auto y0 = model.get_initial_values();
    model.get_flows(y0, y0, 0, dydt_default);

    EXPECT_NEAR(dydt_default[0], 25, 1e-12);
    EXPECT_NEAR(dydt_default[1], 50, 1e-12);
    EXPECT_NEAR(dydt_default[2], 25, 1e-12);
}

TEST(TestSeir, get_flows_two_agegroups)
{
    mio::oseir::Model<double> model(2);
    auto nb_groups = model.parameters.get_num_groups();

    constexpr double total_first_population  = 400;
    constexpr double total_second_population = 200;

    auto& params = model.parameters;

    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations.set_difference_from_group_total<mio::AgeGroup>(
        {mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}, total_first_population);
    model.populations[{mio::AgeGroup(1), mio::oseir::InfectionState::Exposed}]   = 10;
    model.populations[{mio::AgeGroup(1), mio::oseir::InfectionState::Infected}]  = 10;
    model.populations[{mio::AgeGroup(1), mio::oseir::InfectionState::Recovered}] = 10;
    model.populations.set_difference_from_group_total<mio::AgeGroup>(
        {mio::AgeGroup(1), mio::oseir::InfectionState::Susceptible}, total_second_population);

    for (auto i = mio::AgeGroup(0); i <= mio::AgeGroup(1); i++) {
        model.parameters.get<mio::oseir::TimeExposed<double>>()[i]                      = 2;
        model.parameters.get<mio::oseir::TimeInfected<double>>()[i]                     = 4;
        model.parameters.get<mio::oseir::TransmissionProbabilityOnContact<double>>()[i] = 1;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::oseir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 1.0));
    model.check_constraints();

    auto dydt_default = Eigen::VectorXd(6);
    dydt_default.setZero();
    auto y0 = model.get_initial_values();
    model.get_flows(y0, y0, 0, dydt_default);

    EXPECT_NEAR(dydt_default[0], 30, 1e-12);
    EXPECT_NEAR(dydt_default[1], 50, 1e-12);
    EXPECT_NEAR(dydt_default[2], 25, 1e-12);
    EXPECT_NEAR(dydt_default[3], 51, 1e-12);
    EXPECT_NEAR(dydt_default[4], 5, 1e-12);
    EXPECT_NEAR(dydt_default[5], 2.5, 1e-12);
}

TEST(TestSeir, Simulation)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 1;

    mio::oseir::Model<double> model(1);

    constexpr double total_population                                            = 400;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];

    model.parameters.set<mio::oseir::TimeExposed<double>>(2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(4);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(1);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(1);

    model.check_constraints();

    std::shared_ptr<mio::IntegratorCore<double, 1>> integrator = std::make_shared<mio::EulerIntegratorCore<double>>();

    auto sim = simulate(t0, tmax, dt, model, integrator);

    EXPECT_EQ(sim.get_num_time_points(), 2);

    const auto& results_t1 = sim.get_last_value();
    EXPECT_NEAR(results_t1[0], 75, 1e-12);
    EXPECT_NEAR(results_t1[1], 75, 1e-12);
    EXPECT_NEAR(results_t1[2], 125, 1e-12);
    EXPECT_NEAR(results_t1[3], 125, 1e-12);
}

TEST(TestSeir, FlowSimulation)
{
    double t0   = 0;
    double tmax = 1;
    double dt   = 1;

    mio::oseir::Model model(1);

    constexpr double total_population                                            = 400;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}];

    model.parameters.set<mio::oseir::TimeExposed<double>>(2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(4);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(1);
    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(1);

    model.check_constraints();

    std::shared_ptr<mio::IntegratorCore<double, 1>> integrator = std::make_shared<mio::EulerIntegratorCore<double>>();

    auto sim = simulate_flows(t0, tmax, dt, model, integrator);

    // results
    EXPECT_EQ(sim[0].get_num_time_points(), 2);
    const auto& results_t1 = sim[0].get_last_value();
    EXPECT_NEAR(results_t1[0], 75, 1e-12);
    EXPECT_NEAR(results_t1[1], 75, 1e-12);
    EXPECT_NEAR(results_t1[2], 125, 1e-12);
    EXPECT_NEAR(results_t1[3], 125, 1e-12);

    // flows
    EXPECT_EQ(sim[1].get_num_time_points(), 2);
    const auto& flows_t1 = sim[1].get_last_value();
    EXPECT_NEAR(flows_t1[0], 25, 1e-12);
    EXPECT_NEAR(flows_t1[1], 50, 1e-12);
    EXPECT_NEAR(flows_t1[2], 25, 1e-12);
}
