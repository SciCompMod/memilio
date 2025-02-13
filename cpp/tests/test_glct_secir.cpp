/* 
* Copyright (C) 2020-2025 MEmilio
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

#include "glct_secir/model.h"
#include "glct_secir/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/eigen.h"
#include "load_test_data.h"

#include <gtest/gtest.h>
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"

// Test if the function eval_right_hand_side() is working using a hand calculated result.
TEST(TestGLCTSecir, testEvalRightHandSide)
{
    // Define initial values, parameters and numbers of subcompartments according to the choices of the
    // testEvalRightHandSide of the LCT testing suite. For more details,
    // we refer to the example glct_secir.cpp.
    using Model          = mio::glsecir::Model<2, 6, 4, 4, 4>;
    using LctState       = Model::LctState;
    using InfectionState = LctState::InfectionState;

    Model model;

    // Set parameters such that the stay times are Erlang-distributed as in the corresponding LCT model.
    // Exposed.
    // Default functions are used to set the parameters but the corresponding dimensions have to be set manually.
    model.parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>());
    model.parameters.get<mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms>() =
        mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>(), 3.2);
    // InfectedNoSymptoms.
    Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedNoSymptoms = Eigen::VectorX<ScalarType>::Zero(
        (Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>());
    StartingProbabilitiesInfectedNoSymptoms[0]                                         = 1 - 0.09;
    StartingProbabilitiesInfectedNoSymptoms[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.)] = 0.09;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        StartingProbabilitiesInfectedNoSymptoms;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.), 2.);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.), 2.);
    // InfectedSymptoms.
    Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedSymptoms = Eigen::VectorX<ScalarType>::Zero(
        (Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>());
    StartingProbabilitiesInfectedSymptoms[0]                                         = 0.2;
    StartingProbabilitiesInfectedSymptoms[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.)] = 1 - 0.2;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() = StartingProbabilitiesInfectedSymptoms;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.), 5.8);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.), 5.8);
    // InfectedSevere.
    Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedSevere = Eigen::VectorX<ScalarType>::Zero(
        (Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedSevere>());
    StartingProbabilitiesInfectedSevere[0]                                         = 0.25;
    StartingProbabilitiesInfectedSevere[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.)] = 1 - 0.25;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() = StartingProbabilitiesInfectedSevere;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical>() =
        mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.), 9.5);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSevereToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.), 9.5);
    // InfectedCritical.
    Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedCritical = Eigen::VectorX<ScalarType>::Zero(
        (Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedCritical>());
    StartingProbabilitiesInfectedCritical[0]                                         = 0.3;
    StartingProbabilitiesInfectedCritical[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.)] = 1 - 0.3;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() = StartingProbabilitiesInfectedCritical;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToDead>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToDead().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.), 7.1);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.), 7.1);

    model.parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::glsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));

    model.parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model.parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    model.parameters.get<mio::glsecir::Seasonality>()                    = 0.;
    model.parameters.get<mio::glsecir::StartDay>()                       = 0;

    // Define initial population distribution in infection states, one entry per subcompartment.
    std::vector<std::vector<ScalarType>> initial_populations = {
        {750},
        {30, 20},
        {20 * StartingProbabilitiesInfectedNoSymptoms[0], 10 * StartingProbabilitiesInfectedNoSymptoms[0],
         10 * StartingProbabilitiesInfectedNoSymptoms[0], 20 * (1 - StartingProbabilitiesInfectedNoSymptoms[0]),
         10 * (1 - StartingProbabilitiesInfectedNoSymptoms[0]), 10 * (1 - StartingProbabilitiesInfectedNoSymptoms[0])},
        {30 * StartingProbabilitiesInfectedSymptoms[0], 20 * StartingProbabilitiesInfectedSymptoms[0],
         30 * (1 - StartingProbabilitiesInfectedSymptoms[0]), 20 * (1 - StartingProbabilitiesInfectedSymptoms[0])},
        {40 * StartingProbabilitiesInfectedSevere[0], 10 * StartingProbabilitiesInfectedSevere[0],
         40 * (1 - StartingProbabilitiesInfectedSevere[0]), 10 * (1 - StartingProbabilitiesInfectedSevere[0])},
        {10 * StartingProbabilitiesInfectedCritical[0], 20 * StartingProbabilitiesInfectedCritical[0],
         10 * (1 - StartingProbabilitiesInfectedCritical[0]), 20 * (1 - StartingProbabilitiesInfectedCritical[0])},
        {20},
        {10}};
    std::vector<ScalarType> flat_initial_populations;
    for (auto&& vec : initial_populations) {
        flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
    }
    Eigen::VectorX<ScalarType> pop(LctState::Count);
    for (size_t i = 0; i < LctState::Count; i++) {
        pop[i] = flat_initial_populations[i];
    }

    // Compare the result of get_derivatives() with a hand calculated result.
    Eigen::VectorX<ScalarType> dydt(LctState::Count);
    model.get_derivatives(pop, pop, 0, dydt);
    // This vector is the equivalent of the result defined in the test suite testEvalRightHandSide of the LCT model.
    Eigen::VectorX<ScalarType> compare(LctState::Count);
    compare << -15.3409, -3.4091, 6.25, -17.5 * 0.91, 15 * 0.91, 0 * 0.91, -17.5 * 0.09, 15 * 0.09, 0 * 0.09,
        3.3052 * 0.2, 3.4483 * 0.2, 3.3052 * 0.8, 3.4483 * 0.8, -7.0417 * 0.25, 6.3158 * 0.25, -7.0417 * 0.75,
        6.3158 * 0.75, -2.2906 * 0.3, -2.8169 * 0.3, -2.2906 * 0.7, -2.8169 * 0.7, 12.3899, 1.6901;

    for (size_t i = 0; i < LctState::Count; i++) {
        EXPECT_NEAR(compare[i], dydt[i], 1e-3) << "Condition failed at index: " << i;
    }
}

// Model setup to compare result with a previous output of an LCT model.
class ModelTestGLCTSecir : public testing::Test
{
public:
    using Model          = mio::glsecir::Model<2, 6, 2, 2, 10>;
    using LctState       = Model::LctState;
    using InfectionState = LctState::InfectionState;

protected:
    virtual void SetUp()
    {
        model = new Model();
        // --- Set parameters. ---
        // Exposed.
        model->parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
            mio::glsecir::StartingProbabilitiesExposed().get_default(
                LctState::get_num_subcompartments<InfectionState::Exposed>());
        model->parameters.get<mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms>() =
            mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms().get_default(
                LctState::get_num_subcompartments<InfectionState::Exposed>(), 3.2);
        // InfectedNoSymptoms.
        Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedNoSymptoms = Eigen::VectorX<ScalarType>::Zero(
            (Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>());
        StartingProbabilitiesInfectedNoSymptoms[0]                                         = 1 - 0.09;
        StartingProbabilitiesInfectedNoSymptoms[(Eigen::Index)(
            LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.)] = 0.09;
        model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
            StartingProbabilitiesInfectedNoSymptoms;
        model->parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>() =
            mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms().get_default(
                (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.), 2.);
        model->parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered>() =
            mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered().get_default(
                (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.), 2.);
        // InfectedSymptoms.
        Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedSymptoms = Eigen::VectorX<ScalarType>::Zero(
            (Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>());
        StartingProbabilitiesInfectedSymptoms[0]                                         = 0.2;
        StartingProbabilitiesInfectedSymptoms[(Eigen::Index)(
            LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.)] = 1 - 0.2;
        model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() =
            StartingProbabilitiesInfectedSymptoms;
        model->parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere>() =
            mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere().get_default(
                (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.), 5.8);
        model->parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered>() =
            mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(
                (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.), 5.8);
        // InfectedSevere.
        Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedSevere = Eigen::VectorX<ScalarType>::Zero(
            (Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedSevere>());
        StartingProbabilitiesInfectedSevere[0]                                         = 0.25;
        StartingProbabilitiesInfectedSevere[(Eigen::Index)(
            LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.)] = 1 - 0.25;
        model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() =
            StartingProbabilitiesInfectedSevere;
        model->parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical>() =
            mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical().get_default(
                (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.), 9.5);
        model->parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToRecovered>() =
            mio::glsecir::TransitionMatrixInfectedSevereToRecovered().get_default(
                (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.), 9.5);
        // InfectedCritical.
        Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedCritical = Eigen::VectorX<ScalarType>::Zero(
            (Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedCritical>());
        StartingProbabilitiesInfectedCritical[0]                                         = 0.3;
        StartingProbabilitiesInfectedCritical[(Eigen::Index)(
            LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.)] = 1 - 0.3;
        model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() =
            StartingProbabilitiesInfectedCritical;
        model->parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToDead>() =
            mio::glsecir::TransitionMatrixInfectedCriticalToDead().get_default(
                (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.), 7.1);
        model->parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>() =
            mio::glsecir::TransitionMatrixInfectedCriticalToRecovered().get_default(
                (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.), 7.1);

        model->parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() = 0.05;

        mio::ContactMatrixGroup& contact_matrix = model->parameters.get<mio::glsecir::ContactPatterns>();
        contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
        contact_matrix[0].add_damping(0.7, mio::SimulationTime(2.));

        model->parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() = 0.7;
        model->parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
        model->parameters.get<mio::glsecir::Seasonality>()                    = 0.;
        model->parameters.get<mio::glsecir::StartDay>()                       = 0;

        // --- Define initial distribution of the population in the subcompartments. ---
        std::vector<std::vector<ScalarType>> initial_populations = {
            {750},
            {30, 20},
            {20 * StartingProbabilitiesInfectedNoSymptoms[0], 10 * StartingProbabilitiesInfectedNoSymptoms[0],
             10 * StartingProbabilitiesInfectedNoSymptoms[0], 20 * (1 - StartingProbabilitiesInfectedNoSymptoms[0]),
             10 * (1 - StartingProbabilitiesInfectedNoSymptoms[0]),
             10 * (1 - StartingProbabilitiesInfectedNoSymptoms[0])},
            {50 * StartingProbabilitiesInfectedSymptoms[0], 50 * (1 - StartingProbabilitiesInfectedSymptoms[0])},
            {50 * StartingProbabilitiesInfectedSevere[0], 50 * (1 - StartingProbabilitiesInfectedSevere[0])},
            {10 * StartingProbabilitiesInfectedCritical[0], 10 * StartingProbabilitiesInfectedCritical[0],
             5 * StartingProbabilitiesInfectedCritical[0], 3 * StartingProbabilitiesInfectedCritical[0],
             2 * StartingProbabilitiesInfectedCritical[0], 10 * (1 - StartingProbabilitiesInfectedCritical[0]),
             10 * (1 - StartingProbabilitiesInfectedCritical[0]), 5 * (1 - StartingProbabilitiesInfectedCritical[0]),
             3 * (1 - StartingProbabilitiesInfectedCritical[0]), 2 * (1 - StartingProbabilitiesInfectedCritical[0])},
            {20},
            {10}};
        // Transfer the initial values in initial_populations to the model->
        std::vector<ScalarType> flat_initial_populations;
        for (auto&& vec : initial_populations) {
            flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
        }
        for (size_t i = 0; i < LctState::Count; i++) {
            model->populations[mio::Index<LctState>(i)] = flat_initial_populations[i];
        }
    }

    virtual void TearDown()
    {
        delete model;
    }

public:
    Model* model = nullptr;
};

// Test compares a simulation with a previous output of an LCT model stored in a .csv file.
// This tests that the GLCT model reduces to an LCT model with the default parameters as well as nothing has changed
// in the implementation such that we still obtain an previous result.
TEST_F(ModelTestGLCTSecir, compareWithPreviousRun)
{
    ScalarType tmax                    = 3;
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, ModelTestGLCTSecir::Model>(
        0, tmax, 0.5, *model,
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>());

    // Compare InfectionState compartments.
    mio::TimeSeries<ScalarType> population = model->calculate_compartments(result);
    auto compare_population                = load_test_data_csv<ScalarType>("lct-secir-compartments-compare.csv");

    ASSERT_EQ(compare_population.size(), static_cast<size_t>(population.get_num_time_points()));
    for (size_t i = 0; i < compare_population.size(); i++) {
        ASSERT_EQ(compare_population[i].size(), static_cast<size_t>(population.get_num_elements()) + 1)
            << "at row " << i;
        EXPECT_NEAR(population.get_time(i), compare_population[i][0], 1e-3) << "at row " << i;
        for (size_t j = 1; j < compare_population[i].size(); j++) {
            EXPECT_NEAR(population.get_value(i)[j - 1], compare_population[i][j], 1e-3) << " at row " << i;
        }
    }
}

// Check constraints of Parameters class.
TEST_F(ModelTestGLCTSecir, testConstraintsModel)
{
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // --- Check correct setup. ---
    bool constraint_check = model->check_constraints();
    EXPECT_FALSE(constraint_check);

    // Check if the number of subcompartments does not match the dimension of the vector with StartingProbabilities.
    Eigen::VectorX<ScalarType> wrong_size = Eigen::VectorX<ScalarType>::Zero(3);
    wrong_size[0]                         = 1;
    // Exposed.
    model->parameters.get<mio::glsecir::StartingProbabilitiesExposed>() = wrong_size;
    constraint_check                                                    = model->check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>());
    // InfectedNoSymptoms.
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() = wrong_size;
    constraint_check                                                               = model->check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        mio::glsecir::StartingProbabilitiesInfectedNoSymptoms().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>());
    // InfectedSymptoms.
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() = wrong_size;
    constraint_check                                                             = model->check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() =
        mio::glsecir::StartingProbabilitiesInfectedSymptoms().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>());
    // InfectedSevere.
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() = wrong_size;
    constraint_check                                                           = model->check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() =
        mio::glsecir::StartingProbabilitiesInfectedSevere().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedSevere>());
    // InfectedCritical.
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() = wrong_size;
    constraint_check                                                             = model->check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() =
        mio::glsecir::StartingProbabilitiesInfectedCritical().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedCritical>());
}

// Check constraints of Parameters class.
TEST_F(ModelTestGLCTSecir, testConstraintsParameters)
{
    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // --- Check correct setup. ---
    bool constraint_check = model->parameters.check_constraints();
    EXPECT_FALSE(constraint_check);

    // --- Parameters affecting the transmission of the virus. ---
    // Check TransmissionProbabilityOnContact.
    model->parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() = 5.1;
    constraint_check                                                        = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() = 0.05;

    // Check RelativeTransmissionNoSymptoms.
    model->parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() = -0.05;
    constraint_check                                                      = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() = 0.7;

    // Check RiskOfInfectionFromSymptomatic.
    model->parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() = 1.1;
    constraint_check                                                      = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() = 0.25;

    // Check Seasonality.
    model->parameters.get<mio::glsecir::Seasonality>() = 0.6;
    constraint_check                                   = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::Seasonality>() = 0.;

    // --- Check with incorrect dimensions. ---
    // Check non-quadratic TransitionMatrixInfectedCriticalToRecovered.
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>() = Eigen::MatrixXd::Zero(2, 3);
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.), 7.1);

    // Check non matching dimensions of TransitionMatrix and vector with StartingProbabilities.
    Eigen::VectorX<ScalarType> wrong_size = Eigen::VectorX<ScalarType>::Zero(3);
    wrong_size[0]                         = 1;
    // Exposed.
    model->parameters.get<mio::glsecir::StartingProbabilitiesExposed>() = wrong_size;
    constraint_check                                                    = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>());
    // InfectedNoSymptoms.
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() = wrong_size;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        mio::glsecir::StartingProbabilitiesInfectedNoSymptoms().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>());
    // InfectedSymptoms.
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() = wrong_size;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() =
        mio::glsecir::StartingProbabilitiesInfectedSymptoms().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>());
    // InfectedSevere.
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() = wrong_size;
    constraint_check                                                           = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() =
        mio::glsecir::StartingProbabilitiesInfectedSevere().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedSevere>());
    // InfectedCritical.
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() = wrong_size;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() =
        mio::glsecir::StartingProbabilitiesInfectedCritical().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedCritical>());

    // --- Check constraints of the starting probability vectors. ---
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>()[1] = 1.5;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>()[0] = 1.1;
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>()[1] = -0.1;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>()[0] = 1.;
    model->parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>()[1] = 0.;

    // --- Check with invalid transition matrices. ---
    // ExposedToInfectedNoSymptoms.
    model->parameters.get<mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms>()(1, 0) = 10;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms>()(1, 0) = 0.1;
    // InfectedNoSymptomsToInfectedSymptoms.
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>()(2, 1) = 50;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>()(2, 1) = 0.1;
    // InfectedNoSymptomsToRecovered.
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered>()(1, 1) = -1.45;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered>()(1, 1) = -1.5;
    // InfectedSymptomsToInfectedSevere.
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere>()(0, 0) = 1.;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere>()(0, 0) = -1.;
    // InfectedSymptomsToRecovered.
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered>()(0, 0) = 0.1;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered>()(0, 0) = -0.1;
    // InfectedSevereToInfectedCritical.
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical>()(0, 0) = 0.01;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical>()(0, 0) = -0.01;
    // InfectedSevereToRecovered.
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToRecovered>()(0, 0) = 50;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToRecovered>()(0, 0) = -0.1;
    // InfectedCriticalToDead.
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToDead>()(3, 1) = 6;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToDead>()(3, 1) = 0.;
    // InfectedCriticalToRecovered.
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>()(0, 4) = 3;
    constraint_check = model->parameters.check_constraints();
    EXPECT_TRUE(constraint_check);
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>()(0, 0) = -3.;
    model->parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>()(0, 4) = 1.2;

    // --- Check with correct parameters. ---
    constraint_check = model->parameters.check_constraints();
    EXPECT_FALSE(constraint_check);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}

// Test calculate_compartments with a TimeSeries that has an incorrect number of elements.
TEST_F(ModelTestGLCTSecir, testCalculatePopWrongSize)
{
    // Deactivate temporarily log output because an error is expected.
    mio::set_log_level(mio::LogLevel::off);
    // TimeSeries has to have LctState::Count elements.
    size_t wrong_size = LctState::Count - 2;
    // Define TimeSeries with wrong_size elements.
    mio::TimeSeries<ScalarType> wrong_num_elements(wrong_size);
    Eigen::VectorX<ScalarType> vec_wrong_size = Eigen::VectorX<ScalarType>::Ones(wrong_size);
    wrong_num_elements.add_time_point(-10, vec_wrong_size);
    wrong_num_elements.add_time_point(-9, vec_wrong_size);
    // Call the calculate_compartments function with the TimeSeries with a wrong number of elements.
    mio::TimeSeries<ScalarType> population = model->calculate_compartments(wrong_num_elements);
    // A TimeSeries of the right size with values -1 is expected.
    ASSERT_EQ(1, population.get_num_time_points());
    for (int i = 0; i < population.get_num_elements(); i++) {
        EXPECT_EQ(-1, population.get_last_value()[i]);
    }
    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}
