/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Lena Ploetzke, Anna Wendler
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
#include "ide_secir/infection_state.h"
#include "ide_secir/model_ide.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <gtest/gtest.h>

class ModelTestIdeSecir : public testing::Test
{
protected:
    virtual void SetUp()
    {
        using Vec = mio::TimeSeries<ScalarType>::Vector;

        //Set initial conditions
        ScalarType N      = 10000;
        ScalarType deaths = 13.10462213;

        int num_transitions = (int)mio::isecir::InfectionTransition::Count;

        Vec vec_init(num_transitions);
        mio::TimeSeries<ScalarType> init(num_transitions);
        vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 25.0;
        vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 1.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 1.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 1.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 1.0;
        vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 1.0;
        init.add_time_point(-10.0, vec_init);
        while (init.get_last_time() < 0) {
            vec_init *= 1.01;
            init.add_time_point(init.get_last_time() + dt, vec_init);
        }

        // Initialize model
        model = new mio::isecir::Model(std::move(init), N, deaths);

        // Set working parameters.
        mio::SmootherCosine smoothcos(2.0);
        mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
        std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        model->parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);
        mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
        contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
        model->parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        mio::ExponentialSurvivalFunction exponential(0.5);
        mio::StateAgeFunctionWrapper prob(exponential);
        model->parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
        model->parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
        model->parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);
        model->parameters.set<mio::isecir::Seasonality>(0.);
        model->parameters.set<mio::isecir::StartDay>(0);

        model->set_tol_for_support_max(1e-10);
    }

    virtual void TearDown()
    {
        delete model;
    }

public:
    mio::isecir::Model* model = nullptr;
    ScalarType dt             = 1;
};

// Check if population stays constant over course of simulation.
TEST_F(ModelTestIdeSecir, checkPopulationConservation)
{
    mio::TimeSeries<ScalarType> compartments = simulate(15, dt, *model);

    ScalarType num_persons_before = 0.0;
    ScalarType num_persons_after  = 0.0;

    for (auto i = 0; i < compartments[0].size(); i++) {
        num_persons_before += compartments[0][i];
        num_persons_after += compartments.get_last_value()[i];
    }

    EXPECT_NEAR(num_persons_after, num_persons_before, 1e-10);
}

// Compare compartments with previous run.
TEST_F(ModelTestIdeSecir, compareWithPreviousRun)
{
    auto compare                             = load_test_data_csv<ScalarType>("ide-secir-compare.csv");
    mio::TimeSeries<ScalarType> compartments = simulate(5, dt, *model);

    ASSERT_EQ(compare.size(), static_cast<size_t>(compartments.get_num_time_points()));
    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(compartments.get_num_elements()) + 1) << "at row " << i;
        ASSERT_NEAR(compartments.get_time(i), compare[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            ASSERT_NEAR(compartments.get_value(i)[j - 1], compare[i][j], 1e-7) << " at row " << i;
        }
    }
}

// Compare transitions with previous run.
TEST_F(ModelTestIdeSecir, compareWithPreviousRunTransitions)
{
    auto compare = load_test_data_csv<ScalarType>("ide-secir-transitions-compare.csv");

    mio::isecir::Simulation sim(*model, dt);
    sim.advance(5);

    auto transitions = sim.get_transitions();

    size_t iter_0 = 0;
    while (transitions.get_time(iter_0) < compare[0][0]) {
        iter_0++;
    }

    for (size_t i = 0; i < compare.size(); i++) {
        ASSERT_EQ(compare[i].size(), static_cast<size_t>(transitions.get_num_elements()) + 1) << "at row " << i;
        ASSERT_NEAR(transitions.get_time(i + iter_0), compare[i][0], 1e-7) << "at row " << i;
        for (size_t j = 1; j < compare[i].size(); j++) {
            ASSERT_NEAR(transitions.get_value(i + iter_0)[j - 1], compare[i][j], 1e-7) << " at row " << i;
        }
    }
}

// Check that the start time of the simulation is determined by the given time points for the transitions.
TEST(IdeSecir, checkStartTime)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax   = 3.0;
    ScalarType N      = 10000.;
    ScalarType deaths = 10.;
    ScalarType dt     = 1.;
    ScalarType t0     = 2.;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_transitions);

    // Define transitions that will be used for initialization.
    Vec vec_init                                                          = Vec::Constant(num_transitions, 0.);
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed] = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    // Add initial time point to TimeSeries.
    init.add_time_point(0.0, vec_init);
    // Add further time points until t0.
    while (init.get_last_time() < t0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, deaths);

    // Create class for simulation.
    mio::isecir::Simulation sim(model, dt);

    // Check that the last time point of transitions is equal to t0.
    mio::TimeSeries<ScalarType> transitions = sim.get_transitions();
    EXPECT_NEAR(t0, transitions.get_last_time(), 1e-8);

    // Carry out simulation and check that first time point of resulting compartments is equal to t0.
    sim.advance(tmax);
    mio::TimeSeries<ScalarType> compartments = sim.get_result();
    EXPECT_NEAR(t0, compartments.get_time(0), 1e-8);
}

// Check results of our simulation with an example calculated by hand,
// for calculations see internal Overleaf document.
// TODO: Add link to material when published.
TEST(IdeSecir, checkSimulationFunctions)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax   = 0.5;
    ScalarType N      = 10000;
    ScalarType deaths = 10;
    ScalarType dt     = 0.5;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_transitions);

    // Add time points for initialization for transitions and death.
    Vec vec_init                                                          = Vec::Constant(num_transitions, 0.);
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed] = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    // Add initial time point to TimeSeries.
    init.add_time_point(-0.5, vec_init);
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, deaths);

    // Set working parameters.
    // In this example, SmootherCosine with parameter 1 (and thus with a maximum support of 1)
    // is used for all TransitionDistribution%s.
    mio::SmootherCosine smoothcos(1.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.5);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 4.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::SmootherCosine smoothcos_prob(1.0);
    mio::StateAgeFunctionWrapper prob(smoothcos_prob);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);
    model.parameters.set<mio::isecir::Seasonality>(0.);
    model.parameters.set<mio::isecir::StartDay>(0);

    // Carry out simulation.
    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated    = sim.get_result();
    mio::TimeSeries<ScalarType> transitions_simulated = sim.get_transitions();

    // Define vectors for compartments and transitions with values from example
    // (calculated by hand, see internal Overleaf document).
    // TODO: Add link to material when published.
    Vec secihurd0((int)mio::isecir::InfectionState::Count);
    Vec secihurd1((int)mio::isecir::InfectionState::Count);
    Vec transitions1(num_transitions);
    secihurd0 << 4995, 0.5, 0, 4, 0, 0, 4990.5, 10;
    secihurd1 << 4994.00020016, 0.49989992, 0.49994996, 0.12498749, 1.03124687, 0.25781172, 4993.45699802, 10.12890586;
    transitions1 << 0.99979984, 0.99989992, 0.24997498, 0.24997498, 2.06249374, 2.06249374, 0.51562344, 0.51562344,
        0.12890586, 0.12890586;

    // Compare SECIHURD compartments at times 0 and 1.
    for (Eigen::Index i = 0; i < (Eigen::Index)mio::isecir::InfectionState::Count; i++) {
        EXPECT_NEAR(secihurd_simulated[0][i], secihurd0[i], 1e-8);
        EXPECT_NEAR(secihurd_simulated[1][i], secihurd1[i], 1e-8);
    }

    // Compare transitions at time 1.
    for (Eigen::Index i = 0; i < num_transitions; i++) {
        EXPECT_NEAR(transitions_simulated[transitions_simulated.get_num_time_points() - 1][i], transitions1[i], 1e-8);
    }
}

// Check if the model uses the correct method for initialization using the function get_initialization_method_compartments().
TEST(IdeSecir, checkInitializations)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax   = 1;
    ScalarType N      = 10000;
    ScalarType deaths = 13.10462213;
    ScalarType dt     = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_transitions);
    // Add initial time point to time series.
    init.add_time_point(-10, Vec::Constant(num_transitions, 3.0));
    // Add further time points until time 0.
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, Vec::Constant(num_transitions, 3.0));
    }

    // --- Case with total_confirmed_cases.
    mio::TimeSeries<ScalarType> init_copy1(init);
    mio::isecir::Model model1(std::move(init_copy1), N, deaths, 1000);

    // Check that the initialization method is not already set.
    EXPECT_EQ(0, model1.get_initialization_method_compartments());

    // Carry out simulation.
    mio::isecir::Simulation sim1(model1, dt);
    sim1.advance(tmax);

    // Verify that the expected initialization method was used.
    EXPECT_EQ(1, sim1.get_model().get_initialization_method_compartments());

    // --- Case with S.
    /* !! For the other tests, the contact rate is set to 0 so that the force of infection is zero.
     The forceofinfection initialization method is therefore not used for these tests.*/
    mio::isecir::Parameters parameters;
    mio::ContactMatrixGroup contact_matrix         = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                              = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 0));
    parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::TimeSeries<ScalarType> init_copy2(init);
    mio::isecir::Model model2(std::move(init_copy2), N, deaths, 0, std::move(parameters));

    model2.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Susceptible] = 5000;
    model2.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered]   = 0;

    // Carry out simulation.
    mio::isecir::Simulation sim2(model2, dt);
    sim2.advance(tmax);

    // Verify that the expected initialization method was used.
    EXPECT_EQ(2, sim2.get_model().get_initialization_method_compartments());

    // --- Case with R.
    mio::TimeSeries<ScalarType> init_copy3(init);
    mio::isecir::Model model3(std::move(init_copy3), N, deaths, 0, std::move(parameters));

    model3.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Susceptible] = 0;
    model3.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered]   = 1000;

    // Carry out simulation.
    mio::isecir::Simulation sim3(model3, dt);
    sim3.advance(tmax);

    // Verify that the expected initialization method was used.
    EXPECT_EQ(3, sim3.get_model().get_initialization_method_compartments());

    // --- Case with forceofinfection.
    mio::TimeSeries<ScalarType> init_copy4(init);
    mio::isecir::Model model4(std::move(init_copy4), N, deaths, 0);

    // Carry out simulation.
    mio::isecir::Simulation sim4(model4, dt);
    sim4.advance(tmax);

    // Verify that the expected initialization method was used.
    EXPECT_EQ(4, sim4.get_model().get_initialization_method_compartments());

    // --- Case without fitting initialization method.
    // Deactivate temporarily log output for next tests. Errors are expected here.
    mio::set_log_level(mio::LogLevel::off);

    mio::TimeSeries<ScalarType> init_copy5(init);
    mio::isecir::Model model5(std::move(init_copy5), N, deaths, 0, std::move(parameters));

    model5.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Susceptible] = 0;
    model5.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered]   = 0;

    // Carry out simulation.
    mio::isecir::Simulation sim5(model5, dt);
    sim5.advance(tmax);

    // Verify that initialization was not possible with one of the models methods.
    EXPECT_EQ(-1, sim5.get_model().get_initialization_method_compartments());

    // --- Test with negative number of deaths.
    deaths = -10;

    // Here we do not need a copy of init as this is the last use of the vector. We can apply move directly.
    mio::isecir::Model model6(std::move(init), N, deaths, 0, std::move(parameters));

    model6.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Susceptible] = 0;
    model6.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered]   = 0;

    // Carry out simulation.
    mio::isecir::Simulation sim6(model6, dt);
    sim6.advance(tmax);

    // Verify that initialization was possible but the result is not appropriate.
    EXPECT_EQ(-2, sim6.get_model().get_initialization_method_compartments());

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

// a) Test if the function check_constraints() of the class Model correctly reports errors in the model constraints.
// b) Test if check_constraints() does not complain if the conditions are met.
TEST(IdeSecir, testModelConstraints)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Set wrong initial data and use check_constraints().
    // Follow the same order as in check_constraints().

    // --- Test with wrong size of the initial value vector for the flows.
    ScalarType N      = 10000;
    ScalarType deaths = 10;
    ScalarType dt     = 1;

    int num_transitions  = (int)mio::isecir::InfectionTransition::Count;
    int num_compartments = (int)mio::isecir::InfectionState::Count;

    // Create TimeSeries of the wrong size.
    mio::TimeSeries<ScalarType> init_wrong_size(num_transitions + 1);
    // Add time points with vectors of the wrong size.
    Vec vec_init_wrong_size = Vec::Constant(num_transitions + 1, 0.);
    vec_init_wrong_size[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = 10.0;
    init_wrong_size.add_time_point(-3, vec_init_wrong_size);
    while (init_wrong_size.get_last_time() < 0) {
        init_wrong_size.add_time_point(init_wrong_size.get_last_time() + dt, vec_init_wrong_size);
    }

    // Initialize a model.
    mio::isecir::Model model_wrong_size(std::move(init_wrong_size), N, deaths);

    auto constraint_check = model_wrong_size.check_constraints(dt);
    EXPECT_TRUE(constraint_check);

    // --- Test with negative number of deaths.
    // Create TimeSeries with num_transitions elements.
    mio::TimeSeries<ScalarType> init(num_transitions);
    // Add time points for initialization of transitions.
    Vec vec_init                                                                 = Vec::Constant(num_transitions, 0.);
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = 10.0;
    init.add_time_point(-3, vec_init);
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }
    deaths = -10;

    // Initialize a model.
    mio::isecir::Model model_negative_deaths(std::move(init), N, deaths);

    // Return true for negative entry in m_populations.
    constraint_check = model_negative_deaths.check_constraints(dt);
    EXPECT_TRUE(constraint_check);

    // --- Test with too few time points.
    // Create TimeSeries with num_transitions elements.
    mio::TimeSeries<ScalarType> init_few_timepoints(num_transitions);
    // Add time points for initialization of transitions.
    init_few_timepoints.add_time_point(-3, vec_init);
    while (init_few_timepoints.get_last_time() < 0) {
        init_few_timepoints.add_time_point(init_few_timepoints.get_last_time() + dt, vec_init);
    }
    deaths = 10;

    // Initialize a model.
    mio::isecir::Model model(std::move(init_few_timepoints), N, deaths);

    mio::ExponentialSurvivalFunction exponential(4.0);
    mio::StateAgeFunctionWrapper delaydistribution(exponential);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Return true for not enough time points given for the initial transitions.
    constraint_check = model.check_constraints(dt);
    EXPECT_TRUE(constraint_check);

    // --- Test with negative transitions.
    // Create TimeSeries with num_transitions elements.
    mio::TimeSeries<ScalarType> init_negative_transitions(num_transitions);
    // Add time points for initialization of transitions.
    init_negative_transitions.add_time_point(-3, vec_init);
    while (init_negative_transitions.get_last_time() < 0) {
        init_negative_transitions.add_time_point(init_negative_transitions.get_last_time() + dt, (-1) * vec_init);
    }

    // Initialize a model.
    mio::isecir::Model model_negative_transitions(std::move(init_negative_transitions), N, deaths);

    // Return true for negative entries in the initial transitions.
    constraint_check = model_negative_transitions.check_constraints(dt);
    EXPECT_TRUE(constraint_check);

    // --- Test with last time point of transitions not matching last time point of populations.
    // Create TimeSeries with num_transitions elements.
    mio::TimeSeries<ScalarType> init_different_last_time(num_transitions);
    // Add enough time points for initialization of transitions but with different last time point
    // than before so that it does not match last time point of m_populations (that was set in
    // when constructing model above).
    init_different_last_time.add_time_point(-4, vec_init);
    while (init_different_last_time.get_last_time() < 1) {
        init_different_last_time.add_time_point(init_different_last_time.get_last_time() + dt, vec_init);
    }

    model.m_transitions = init_different_last_time;

    // Return true for not last time points of compartments and transitions not matching.
    constraint_check = model.check_constraints(dt);
    EXPECT_TRUE(constraint_check);

    // --- Test with TimeSeries for populations that contains more than one time point.
    // Create TimeSeries with num_compartments elements.
    mio::TimeSeries<ScalarType> populations_many_timepoints(num_compartments);
    // Add time points.
    Vec vec_populations = Vec::Constant(num_compartments, 0.);
    populations_many_timepoints.add_time_point(0, vec_populations);
    while (populations_many_timepoints.get_last_time() < 1) {
        populations_many_timepoints.add_time_point(populations_many_timepoints.get_last_time() + dt, vec_populations);
    }

    model.m_populations = populations_many_timepoints;

    // Return true for too many time points given for populations.
    constraint_check = model.check_constraints(dt);
    EXPECT_TRUE(constraint_check);

    // --- Correct wrong setup so that next check can go through.
    // Create TimeSeries with num_compartments elements.
    mio::TimeSeries<ScalarType> correct_populations(num_compartments);
    // Add one time point.
    correct_populations.add_time_point(1, vec_populations);

    model.m_populations = correct_populations;

    constraint_check = model.check_constraints(dt);
    EXPECT_FALSE(constraint_check);

    // --- The check_constraints() function of parameters is tested in its own test below. ---

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

// a) Test if check_constraints() function of Parameters correctly reports wrongly set parameters.
// b) Test if check_constraints() function of Parameters does not complain if parameters are set within correct ranges.
TEST(IdeSecir, testParametersConstraints)
{
    // Create an object from the class Parameters.
    mio::isecir::Parameters parameters;

    // Deactivate temporarily log output for next tests as warnings are expected.
    mio::set_log_level(mio::LogLevel::off);

    // Set wrong parameters and test if check_constraints() reports them correctly.
    // Test in the same order as in check_constraints().
    // Create invalid and valid function.
    mio::ConstantFunction constant_func_neg(-1);
    mio::StateAgeFunctionWrapper prob_neg(constant_func_neg);
    mio::ConstantFunction constant_func_pos(1);
    mio::StateAgeFunctionWrapper prob_pos(constant_func_pos);

    // Warn, i.e., return true for wrong TransmissionProbabilityOnContact.
    parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob_neg);
    auto constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Correct wrong parameter so that next check can go through.
    parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob_pos);
    // Warn, i.e., return true for wrong RelativeTransmissionNoSymptoms.
    parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob_neg);
    constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Correct wrong parameter so that next check can go through.
    parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob_pos);
    // Warn, i.e., return true for wrong RiskOfInfectionFromSymptomatic.
    parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob_neg);
    constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Correct wrong parameter so that next check can go through.
    parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob_pos);

    // Set wrong values for InfectionTransitions.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]          = 0.2;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)]   = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)]     = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)]   = 0.4;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)]        = -0.6;
    parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    // Negative probability for InfectedCriticalToDead.
    constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check if SusceptibleToExposed is not equal to 1.
    parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        0.6;
    constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check if ExposedToInfectedNoSymptoms is not equal to 1.
    parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::SusceptibleToExposed] = 1.0;
    parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]                                        = 0.6;
    constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check sum InfectedNoSymptomsToInfectedSymptoms + InfectedNoSymptomsToRecovered.
    parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]   = 1.0;
    parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] = 0.9;
    constraint_check                                                         = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check sum InfectedSymptomsToInfectedSevere + InfectedSymptomsToRecovered.
    parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]    = 0.0;
    parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] = 0.2;
    constraint_check                                                            = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check sum InfectedSevereToInfectedCritical + InfectedCriticalToRecovered.
    parameters.get<mio::isecir::TransitionProbabilities>()[(
        int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] = 1.0;
    parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
        0.4;
    constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check sum InfectedCriticalToDead + InfectedSevereToRecovered.
    parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
        0.0;
    parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        0.59;
    constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Set wrong function type with unlimited support.
    parameters
        .get<mio::isecir::TransitionProbabilities>()[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
        0.6;
    mio::ConstantFunction const_func(1.0);
    mio::StateAgeFunctionWrapper delaydistribution(const_func);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib((int)mio::isecir::InfectionTransition::Count,
                                                               delaydistribution);
    parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);
    constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check wrong value for Seasonality.
    mio::ExponentialSurvivalFunction exponential(4.0);
    mio::StateAgeFunctionWrapper delaydistribution2(exponential);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib2((int)mio::isecir::InfectionTransition::Count,
                                                                delaydistribution2);
    parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib2);
    parameters.set<mio::isecir::Seasonality>(2.);
    constraint_check = parameters.check_constraints();
    EXPECT_TRUE(constraint_check);

    // Check if all parameter are correct.
    parameters.set<mio::isecir::Seasonality>(0.1);
    constraint_check = parameters.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}
// The idea of this test is to check whether the proportion between Recovered and Dead is as expected
// (after simulation for a long enough time, i.e. when the equlibrium is reached).
TEST(IdeSecir, checkProportionRecoveredDeath)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax   = 30;
    ScalarType N      = 10000;
    ScalarType deaths = 10;
    ScalarType dt     = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_transitions);

    // Add time points for initialization for transitions.
    Vec vec_init                                                                 = Vec::Constant(num_transitions, 0.);
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] = 10.0;
    // Add initial time point to TimeSeries.
    init.add_time_point(-12, vec_init);
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, deaths);

    // Set working parameters.
    // All TransitionDistribution%s are ExponentialSurvivalFunction functions.
    // For all TransitionDistribution%s init_parameter=2 is used except for InfectedCriticalToRecovered
    // where init_parameter=3 is used.
    mio::ExponentialSurvivalFunction exponential(4.0);
    mio::StateAgeFunctionWrapper delaydistribution(exponential);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_distribution_parameter(
        3.0);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Set probabilities so that all individuals go from Susceptible to InfectedCritical with probability 1,
    // from there they move to Recovered or Dead with probability 0.4 and 0.6, respectively.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)]   = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)]     = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)]   = 0.4;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)]        = 0.6;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 1.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ExponentialSurvivalFunction exponential2(0.5);
    mio::StateAgeFunctionWrapper prob(exponential2);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);
    // Carry out simulation.
    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated = sim.get_result();

    // Check whether equilibrium has been reached, only then the right proportion
    // between Recovered and Dead is expected.
    EXPECT_TRUE(secihurd_simulated[Eigen::Index(tmax / dt - 1)] == secihurd_simulated[Eigen::Index(tmax / dt - 2)]);

    // Check if the compartments E, C, I, H and U are almost empty at the equilibrium.
    EXPECT_NEAR(secihurd_simulated.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Exposed], 0, 1e-18);
    EXPECT_NEAR(secihurd_simulated.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::InfectedNoSymptoms], 0,
                1e-8);
    EXPECT_NEAR(secihurd_simulated.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::InfectedSymptoms], 0,
                1e-8);
    EXPECT_NEAR(secihurd_simulated.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::InfectedSevere], 0,
                1e-8);
    EXPECT_NEAR(secihurd_simulated.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::InfectedCritical], 0,
                1e-8);
                
    // Check whether equilibrium has the right proportion between Recovered and Dead.
    EXPECT_NEAR((vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] /
                 vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)]) *
                    (secihurd_simulated.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered] -
                     secihurd_simulated[0][(Eigen::Index)mio::isecir::InfectionState::Recovered]),
                secihurd_simulated.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Dead] -
                    secihurd_simulated[0][(Eigen::Index)mio::isecir::InfectionState::Dead],
                1e-8);
}

// The idea of this test is to confirm that the equilibrium of the compartments
// (after simulation for a long enough time) does not change if a different TransitionDistribution is used
// for the transition from InfectedCritical To Recovered.
// It is also checked whether the equilibirum is reached at an earlier time if m_support_max is chosen smaller.
TEST(IdeSecir, compareEquilibria)
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax   = 20;
    ScalarType N      = 10000;
    ScalarType deaths = 10;
    ScalarType dt     = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_transitions);

    // Add time points for initialization for transitions.
    Vec vec_init                                                                 = Vec::Constant(num_transitions, 0.);
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = 10.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] = 10.0;
    // Add initial time point to TimeSeries.
    init.add_time_point(-12, vec_init);
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    mio::TimeSeries<ScalarType> init2(init);

    // Initialize two models.
    mio::isecir::Model model(std::move(init), N, deaths);
    mio::isecir::Model model2(std::move(init2), N, deaths);

    // Set working parameters.
    // Here the maximum support for the TransitionDistribution%s is set differently for each model
    // In both models, all TransitionDistribution%s are SmootherCosine.

    // For the Model model.
    // All TransitionDistribution%s have parameter=2.
    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // For the Model model2.
    // All TransitionDistribution%s have parameter=2 except for InfectedCriticalToRecovered
    // which has parameter=7.
    mio::SmootherCosine smoothcos2(2.0);
    mio::StateAgeFunctionWrapper delaydistribution2(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib2(num_transitions, delaydistribution2);
    vec_delaydistrib2[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_distribution_parameter(
        7.0);
    model2.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib2);

    // All remaining parameters are equal for both models.
    // Set probabilities so that all individuals go from Susceptible to InfectedCritical with probability 1,
    // from there they move to Recovered or Dead with probability 0.5, respectively.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)]   = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)]     = 0.0;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)]   = 0.5;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)]        = 0.5;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);
    model2.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix                = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                     = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 1.));
    model.parameters.get<mio::isecir::ContactPatterns>()  = mio::UncertainContactMatrix(contact_matrix);
    model2.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ExponentialSurvivalFunction exponential(0.5);
    mio::StateAgeFunctionWrapper prob(exponential);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);

    model2.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model2.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model2.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);

    // Carry out simulation.
    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated = sim.get_result();

    mio::isecir::Simulation sim2(model2, dt);
    sim2.advance(tmax);
    mio::TimeSeries<ScalarType> secihurd_simulated2 = sim2.get_result();

    // Check whether equilibrium has been reached, only then it makes sense to compare results and times when
    // equilibrium was reached.
    EXPECT_TRUE(secihurd_simulated[Eigen::Index(tmax / dt - 1)] == secihurd_simulated[Eigen::Index(tmax / dt - 2)]);
    EXPECT_TRUE(secihurd_simulated2[Eigen::Index(tmax / dt - 1)] == secihurd_simulated2[Eigen::Index(tmax / dt - 2)]);

    // Check whether both models have the same result at time tmax.
    for (Eigen::Index i = 0; i < (Eigen::Index)mio::isecir::InfectionState::Count; i++) {
        EXPECT_NEAR(secihurd_simulated.get_last_value()[i], secihurd_simulated2.get_last_value()[i], 1e-8);
    }

    // Compute at what time the equilibrium was reached and check whether that time point is smaller for model than
    //for model2 (as a smaller maximum support is used in model compared to model2).
    ScalarType equilibrium_time{};
    ScalarType equilibrium_time2{};
    for (int t = 0; t < secihurd_simulated.get_num_time_points() - 1; t++) {
        if (secihurd_simulated[t] == secihurd_simulated[t + 1]) {
            equilibrium_time = t;
            break;
        }
    }
    for (int t = 0; t < secihurd_simulated.get_num_time_points() - 1; t++) {
        if (secihurd_simulated2[t] == secihurd_simulated2[t + 1]) {
            equilibrium_time2 = t;
            break;
        }
    }

    EXPECT_TRUE(equilibrium_time <= equilibrium_time2);
}

TEST(IdeSecir, checkInfectionTransitions)
{
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.size(), mio::isecir::InfectionTransitionsCount);

    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(0),
              std::make_pair(mio::isecir::InfectionState::Susceptible, mio::isecir::InfectionState::Exposed));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(1),
              std::make_pair(mio::isecir::InfectionState::Exposed, mio::isecir::InfectionState::InfectedNoSymptoms));
    EXPECT_EQ(
        mio::isecir::InfectionTransitionsMap.at(2),
        std::make_pair(mio::isecir::InfectionState::InfectedNoSymptoms, mio::isecir::InfectionState::InfectedSymptoms));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(3),
              std::make_pair(mio::isecir::InfectionState::InfectedNoSymptoms, mio::isecir::InfectionState::Recovered));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(4), std::make_pair(mio::isecir::InfectionState::InfectedSymptoms,
                                                                         mio::isecir::InfectionState::InfectedSevere));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(5),
              std::make_pair(mio::isecir::InfectionState::InfectedSymptoms, mio::isecir::InfectionState::Recovered));
    EXPECT_EQ(
        mio::isecir::InfectionTransitionsMap.at(6),
        std::make_pair(mio::isecir::InfectionState::InfectedSevere, mio::isecir::InfectionState::InfectedCritical));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(7),
              std::make_pair(mio::isecir::InfectionState::InfectedSevere, mio::isecir::InfectionState::Recovered));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(8),
              std::make_pair(mio::isecir::InfectionState::InfectedCritical, mio::isecir::InfectionState::Dead));
    EXPECT_EQ(mio::isecir::InfectionTransitionsMap.at(9),
              std::make_pair(mio::isecir::InfectionState::InfectedCritical, mio::isecir::InfectionState::Recovered));
}
