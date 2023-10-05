/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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

#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/simulation.h"
#include "lct_secir/parameters.h"
#include "ode_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"
#include "load_test_data.h"

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"

// Test confirms that default construction works
TEST(TestLCTSecir, simulateDefault)
{
    ScalarType t0   = 0;
    ScalarType tmax = 1;
    ScalarType dt   = 0.1;

    Eigen::VectorXd init = Eigen::VectorXd::Constant((int)mio::lsecir::InfectionStateBase::Count, 15);
    init[0]              = 200;
    init[3]              = 50;
    init[5]              = 30;

    mio::lsecir::Model model(init);
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(t0, tmax, dt, model);

    EXPECT_NEAR(result.get_last_time(), tmax, 1e-10);
    ScalarType sum_pop = init.sum();
    for (Eigen::Index i = 0; i < result.get_num_time_points(); i++) {
        ASSERT_NEAR(sum_pop, result[i].sum(), 1e-5);
    }
}

/* Test compares the result for an LCT SECIR model with one single Subcompartment for each infection state 
    with the result of the equivalent ODE SECIR model. */
TEST(TestLCTSecir, compareWithOdeSecir)
{
    ScalarType t0   = 0;
    ScalarType tmax = 5;
    ScalarType dt   = 0.1;

    Eigen::VectorXd init = Eigen::VectorXd::Constant((int)mio::lsecir::InfectionStateBase::Count, 15);
    init[0]              = 200;
    init[3]              = 50;
    init[5]              = 30;

    mio::lsecir::Model model_lct(init);
    // Set Parameters of the model
    model_lct.parameters.get<mio::lsecir::TimeExposed>()            = 2 * 4.2 - 5.2;
    model_lct.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2 * (5.2 - 4.2);
    model_lct.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
    model_lct.parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
    model_lct.parameters.get<mio::lsecir::TimeInfectedCritical>()   = 7.1;

    model_lct.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix_lct = model_lct.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix_lct[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    contact_matrix_lct[0].add_damping(0.7, mio::SimulationTime(2.));

    model_lct.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model_lct.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    model_lct.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    model_lct.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
    model_lct.parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
    model_lct.parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.3;

    mio::TimeSeries<ScalarType> result_lct = mio::lsecir::simulate(
        t0, tmax, dt, model_lct,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Initialize ODE model with one single age group
    mio::osecir::Model model_ode(1);

    // Set population
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::Exposed)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::InfectedNoSymptoms)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::InfectedSymptoms)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::InfectedSevere)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::InfectedCritical)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::Recovered)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}] =
        init[Eigen::Index(mio::lsecir::InfectionStateBase::Dead)];
    model_ode.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                    init.sum());

    // set parameters fitting to these of the lct model
    // no restrictions by additional parameters
    model_ode.parameters.set<mio::osecir::StartDay>(0);
    model_ode.parameters.set<mio::osecir::Seasonality>(0);
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity>() = std::numeric_limits<double>::max();
    model_ode.parameters.get<mio::osecir::ICUCapacity>()          = std::numeric_limits<double>::max();

    model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] =
        5.2; // TimeExposed = 2 * SerialInterval - IncubationTime
    model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] =
        4.2; // TimeInfectedNoSymptoms = 2* (IncubationTime - SerialInterval)
    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0] = 5.8;
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]   = 9.5;
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0] = 7.1;

    mio::ContactMatrixGroup& contact_matrix_ode = model_ode.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix_ode[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    contact_matrix_ode[0].add_damping(0.7, mio::SimulationTime(2.));

    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 0.05;
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 0.7;
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]   = 0.09;
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 0.25;
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]        = 0.2;
    model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]                = 0.25;
    model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]                = 0.3;

    mio::TimeSeries<double> result_ode =
        simulate(t0, tmax, dt, model_ode,
                 std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>());

    ASSERT_EQ(result_lct.get_num_time_points(), result_ode.get_num_time_points());
    for (int i = 0; i < 4; ++i) {
        ASSERT_NEAR(result_lct.get_time(i), result_ode.get_time(i), 1e-5);

        ASSERT_NEAR(result_lct[i][(int)mio::lsecir::InfectionStateBase::Susceptible],
                    result_ode[i][(int)mio::osecir::InfectionState::Susceptible], 1e-5);
        ASSERT_NEAR(result_lct[i][(int)mio::lsecir::InfectionStateBase::Exposed],
                    result_ode[i][(int)mio::osecir::InfectionState::Exposed], 1e-5);
        ASSERT_NEAR(result_lct[i][(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms],
                    result_ode[i][(int)mio::osecir::InfectionState::InfectedNoSymptoms], 1e-5);
        ASSERT_NEAR(0, result_ode[i][(int)mio::osecir::InfectionState::InfectedNoSymptomsConfirmed], 1e-5);
        ASSERT_NEAR(result_lct[i][(int)mio::lsecir::InfectionStateBase::InfectedSymptoms],
                    result_ode[i][(int)mio::osecir::InfectionState::InfectedSymptoms], 1e-5);
        ASSERT_NEAR(0, result_ode[i][(int)mio::osecir::InfectionState::InfectedSymptomsConfirmed], 1e-5);
        ASSERT_NEAR(result_lct[i][(int)mio::lsecir::InfectionStateBase::InfectedCritical],
                    result_ode[i][(int)mio::osecir::InfectionState::InfectedCritical], 1e-5);
        ASSERT_NEAR(result_lct[i][(int)mio::lsecir::InfectionStateBase::InfectedSevere],
                    result_ode[i][(int)mio::osecir::InfectionState::InfectedSevere], 1e-5);
        ASSERT_NEAR(result_lct[i][(int)mio::lsecir::InfectionStateBase::Recovered],
                    result_ode[i][(int)mio::osecir::InfectionState::Recovered], 1e-5);
        ASSERT_NEAR(result_lct[i][(int)mio::lsecir::InfectionStateBase::Dead],
                    result_ode[i][(int)mio::osecir::InfectionState::Dead], 1e-5);
    }
}

// test if the right hand side of the model is working
TEST(TestLCTSecir, testEvalRightHandSide)
{ // Setup model

    /* Number of subcompartments, chose more than one subcompartment for all compartments except S, R, D
    so that the function is correct for all selections */
    std::vector<int> SubcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Exposed]            = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 3;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 2;
    mio::lsecir::InfectionState InfState(SubcompartmentNumbers);

    // define initial population distribution in infection states, one entry per subcompartment
    Eigen::VectorXd init(InfState.get_count());
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Susceptible)]            = 750;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Exposed)]                = 30;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Exposed) + 1]            = 20;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms)]     = 20;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 1] = 10;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 2] = 10;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSymptoms)]       = 30;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSymptoms) + 1]   = 20;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSevere)]         = 40;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSevere) + 1]     = 10;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical)]       = 10;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 1]   = 20;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Recovered)]              = 20;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Dead)]                   = 10;

    // initialize model
    mio::lsecir::Model model(std::move(init), InfState);

    // Set parameters of the model
    model.parameters.set<mio::lsecir::TimeExposed>(2 * 4.2 - 5.2);
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2 * (5.2 - 4.2);
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()   = 7.1;

    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    model.parameters.get<mio::lsecir::Seasonality>()                    = 0.1;
    model.parameters.get<mio::lsecir::StartDay>()                       = 0;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.3;

    // compare the result of eval_right_hand_side with a hand calculated result
    size_t num_subcompartments = model.infectionState.get_count();
    Eigen::VectorXd dydt(num_subcompartments);
    model.eval_right_hand_side(model.get_initial_values(), 0, dydt);

    Eigen::VectorXd compare(num_subcompartments);
    compare << -16.875, -1.875, 6.25, -17.5, 15, 0, 3.3052, 3.4483, -7.0417, 6.3158, -2.2906, -2.8169, 12.3899, 1.6901;

    for (size_t i = 0; i < num_subcompartments; i++) {
        ASSERT_NEAR(compare[i], dydt[i], 1e-3);
    }
}

// SetUp for comparing with previous run
class ModelTestLCTSecir : public testing::Test
{
protected:
    virtual void SetUp()
    {
        // Number of subcompartments
        std::vector<int> SubcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
        SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Exposed]            = 2;
        SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 3;
        SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 5;
        mio::lsecir::InfectionState InfState(SubcompartmentNumbers);

        // define initial population distribution in infection states, one entry per Subcompartment
        Eigen::VectorXd init(InfState.get_count());
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Susceptible)]            = 750;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Exposed)]                = 30;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Exposed) + 1]            = 20;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms)]     = 20;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 1] = 10;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 2] = 10;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSymptoms)]       = 50;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSevere)]         = 50;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical)]       = 10;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 1]   = 10;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 2]   = 5;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 3]   = 3;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 4]   = 2;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Recovered)]              = 20;
        init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Dead)]                   = 10;

        // initialize model
        model = new mio::lsecir::Model(std::move(init), InfState);

        // Set Parameters of the model
        model->parameters.get<mio::lsecir::TimeExposed>()            = 2 * 4.2 - 5.2;
        model->parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2 * (5.2 - 4.2);
        model->parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
        model->parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
        model->parameters.get<mio::lsecir::TimeInfectedCritical>()   = 7.1;

        model->parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

        mio::ContactMatrixGroup& contact_matrix = model->parameters.get<mio::lsecir::ContactPatterns>();
        contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
        contact_matrix[0].add_damping(0.7, mio::SimulationTime(2.));

        model->parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
        model->parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
        model->parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
        model->parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
        model->parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
        model->parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.3;
    }

    virtual void TearDown()
    {
        delete model;
    }

public:
    mio::lsecir::Model* model = nullptr;
};

// Test compares a simulation with the result of a previous run stored in a .csv file
TEST_F(ModelTestLCTSecir, compareWithPreviousRun)
{
    ScalarType tmax                    = 3;
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(
        0, tmax, 0.5, *model,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // compare result of subcompartments
    auto compare = load_test_data_csv<ScalarType>("lct-secir-subcompartments-compare.csv");

    ASSERT_EQ(compare.size(), static_cast<size_t>(result.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(result.get_num_elements()) + 1) << "at row " << i;
        ASSERT_NEAR(result.get_time(i), compare[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            ASSERT_NEAR(result.get_value(i)[j - 1], compare[i][j], 1e-7) << " at row " << i;
        }
    }

    // compare results of base compartments
    mio::TimeSeries<ScalarType> population = model->calculate_populations(result);
    auto compare_population                = load_test_data_csv<ScalarType>("lct-secir-compartments-compare.csv");

    ASSERT_EQ(compare_population.size(), static_cast<size_t>(population.get_num_time_points()));
    for (size_t i = 0; i < compare_population.size(); i++) {
        ASSERT_EQ(compare_population[i].size(), static_cast<size_t>(population.get_num_elements()) + 1)
            << "at row " << i;
        ASSERT_NEAR(population.get_time(i), compare_population[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare_population[i].size(); j++) {
            ASSERT_NEAR(population.get_value(i)[j - 1], compare_population[i][j], 1e-7) << " at row " << i;
        }
    }
}

// Test some functions of the model.
TEST_F(ModelTestLCTSecir, testModelFunctions)
{
    mio::TimeSeries<ScalarType> init_wrong_size((int)mio::lsecir::InfectionStateBase::Count);
    Eigen::VectorXd vec_wrong_size = Eigen::VectorXd::Ones((int)mio::lsecir::InfectionStateBase::Count);
    init_wrong_size.add_time_point(-10, vec_wrong_size);
    init_wrong_size.add_time_point(-9, vec_wrong_size);
    mio::TimeSeries<ScalarType> population = model->calculate_populations(init_wrong_size);
    EXPECT_EQ(1, population.get_num_time_points());
    for (int i = 0; i < population.get_num_elements(); i++) {
        EXPECT_EQ(-1, population.get_last_value()[i]);
    }

    EXPECT_TRUE(model->get_heading_CompartmentsBase().compare("S | E | C | I | H | U | R | D") == 0);
    EXPECT_TRUE(model->get_heading_Subcompartments().compare(
                    "S | E1 | E2 | C1 | C2 | C3 | I | H | U1 | U2 | U3 | U4 | U5 | R | D") == 0);

    std::vector<int> SubcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSevere]   = 2;
    mio::lsecir::InfectionState InfState(SubcompartmentNumbers);

    mio::lsecir::Model model2(std::move(Eigen::VectorXd::Ones(InfState.get_count())), InfState);
    EXPECT_TRUE(model2.get_heading_Subcompartments().compare("S | E | C | I1 | I2 | H1 | H2 | U | R | D") == 0);
}

// Check constraints of InfectionState and Parameters.
TEST(TestLCTSecir, testConstraints)
{
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Number of subcompartments
    std::vector<int> SubcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Susceptible]        = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Exposed]            = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 3;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms]   = 4;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedSevere]     = 5;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 6;
    std::vector<int> SubcompartmentNumbers1((int)mio::lsecir::InfectionStateBase::Count - 1, 1);
    mio::lsecir::InfectionState InfState(SubcompartmentNumbers1);

    // Check wrong size of vector.
    bool constraint_check = InfState.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check wrong number of Subcompartments for Susceptibles.
    InfState.set_SubcompartmentNumbers(SubcompartmentNumbers);
    constraint_check = InfState.check_constraints();
    EXPECT_TRUE(constraint_check);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Susceptible] = 1;

    // For Recovered.
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Recovered] = 5;
    InfState.set_SubcompartmentNumbers(SubcompartmentNumbers);
    constraint_check = InfState.check_constraints();
    EXPECT_TRUE(constraint_check);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Recovered] = 1;

    // For Dead.
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Dead] = 3;
    InfState.set_SubcompartmentNumbers(SubcompartmentNumbers);
    constraint_check = InfState.check_constraints();
    EXPECT_TRUE(constraint_check);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Dead] = 1;

    // Check number of Subcompartments is zero.
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Exposed] = 0;
    InfState.set_SubcompartmentNumbers(SubcompartmentNumbers);
    constraint_check = InfState.check_constraints();
    EXPECT_TRUE(constraint_check);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Exposed] = 2;

    // Check with correct parameters.
    InfState.set_SubcompartmentNumbers(SubcompartmentNumbers);
    constraint_check = InfState.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Check for parameters.
    mio::lsecir::Parameters parameters_lct;
    parameters_lct.get<mio::lsecir::TimeExposed>()                      = 0;
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms>()           = 3.1;
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms>()             = 6.1;
    parameters_lct.get<mio::lsecir::TimeInfectedSevere>()               = 11.1;
    parameters_lct.get<mio::lsecir::TimeInfectedCritical>()             = 17.1;
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.01;
    mio::ContactMatrixGroup contact_matrix                              = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0].finalize();
    parameters_lct.get<mio::lsecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 1;
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 1;
    parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.1;
    parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.1;
    parameters_lct.get<mio::lsecir::CriticalPerSevere>()              = 0.1;
    parameters_lct.get<mio::lsecir::DeathsPerCritical>()              = 0.1;

    // Check TimeExposed.
    constraint_check = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeExposed>() = 3.1;

    // Check TimeInfectedNoSymptoms.
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms>() = 0.1;
    constraint_check                                          = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeInfectedNoSymptoms>() = 3.1;

    // Check TimeInfectedSymptoms.
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms>() = -0.1;
    constraint_check                                        = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeInfectedSymptoms>() = 6.1;

    // Check TimeInfectedSevere.
    parameters_lct.get<mio::lsecir::TimeInfectedSevere>() = 0.5;
    constraint_check                                      = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeInfectedSevere>() = 11.1;

    // Check TimeInfectedCritical.
    parameters_lct.get<mio::lsecir::TimeInfectedCritical>() = 0.;
    constraint_check                                        = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TimeInfectedCritical>() = 17.1;

    // Check TransmissionProbabilityOnContact.
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact>() = -1;
    constraint_check                                                    = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.01;

    // Check RelativeTransmissionNoSymptoms.
    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 1.5;
    constraint_check                                                  = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 1;

    // Check RiskOfInfectionFromSymptomatic.
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 1.5;
    constraint_check                                                  = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 1;

    // Check RecoveredPerInfectedNoSymptoms.
    parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 1.5;
    constraint_check                                                  = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.1;

    // Check SeverePerInfectedSymptoms.
    parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>() = -1;
    constraint_check                                             = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::SeverePerInfectedSymptoms>() = 0.1;

    // Check CriticalPerSevere.
    parameters_lct.get<mio::lsecir::CriticalPerSevere>() = -1;
    constraint_check                                     = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::CriticalPerSevere>() = 0.1;

    // Check DeathsPerCritical.
    parameters_lct.get<mio::lsecir::DeathsPerCritical>() = -1;
    constraint_check                                     = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::DeathsPerCritical>() = 0.1;

    // Check Seasonality.
    parameters_lct.get<mio::lsecir::Seasonality>() = 1;
    constraint_check                               = parameters_lct.check_constraints();
    EXPECT_TRUE(constraint_check);
    parameters_lct.get<mio::lsecir::Seasonality>() = 0.1;

    // Check with correct parameters.
    constraint_check = parameters_lct.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Check for model.
    // Check wrong size of initial value vector.
    mio::lsecir::Model model1(std::move(Eigen::VectorXd::Ones(InfState.get_count() - 1)), InfState,
                              std::move(parameters_lct));
    constraint_check = model1.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check with values smaller than zero.
    mio::lsecir::Model model2(std::move(Eigen::VectorXd::Constant(InfState.get_count(), -1)), InfState,
                              std::move(parameters_lct));
    constraint_check = model2.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check with correct conditions.
    mio::lsecir::Model model3(std::move(Eigen::VectorXd::Constant(InfState.get_count(), 100)), InfState,
                              std::move(parameters_lct));
    constraint_check = model3.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}