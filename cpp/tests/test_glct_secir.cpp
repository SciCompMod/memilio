/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "glct_secir/infection_state.h"
#include "glct_secir/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/compartments/simulation.h"
#include "memilio/math/eigen.h"

#include <gtest/gtest.h>

// Test if the function eval_right_hand_side() is working using a hand calculated result.
TEST(TestGLCTSecir, testEvalRightHandSide)
{
    // Define model.
    // Chose more than one subcompartment for all compartments except S, R, D so that the function is correct for all selections.
    using Model          = mio::glsecir::Model<2, 6, 4, 4, 4>;
    using LctState       = Model::LctState;
    using InfectionState = LctState::InfectionState;

    Model model;

    // Set parameters.
    // Exposed.
    model.parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>());
    model.parameters.get<mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms>() =
        mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>(), 3.2);
    // InfectedNoSymptoms.
    Eigen::VectorXd StartingProbabilitiesInfectedNoSymptoms =
        Eigen::VectorXd::Zero((Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>());
    StartingProbabilitiesInfectedNoSymptoms[0] = 1 - 0.09;
    StartingProbabilitiesInfectedNoSymptoms
        [(Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.] = 0.09;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        StartingProbabilitiesInfectedNoSymptoms;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2., 2.);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2., 2.);
    // InfectedSymptoms.
    Eigen::VectorXd StartingProbabilitiesInfectedSymptoms =
        Eigen::VectorXd::Zero((Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>());
    StartingProbabilitiesInfectedSymptoms[0] = 0.2;
    StartingProbabilitiesInfectedSymptoms
        [(Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.] = 1 - 0.2;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() = StartingProbabilitiesInfectedSymptoms;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2., 5.8);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2., 5.8);
    // InfectedSevere.
    Eigen::VectorXd StartingProbabilitiesInfectedSevere =
        Eigen::VectorXd::Zero((Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedSevere>());
    StartingProbabilitiesInfectedSevere[0] = 0.25;
    StartingProbabilitiesInfectedSevere
        [(Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.] = 1 - 0.25;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() = StartingProbabilitiesInfectedSevere;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical>() =
        mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2., 9.5);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSevereToRecovered().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2., 9.5);
    // InfectedCritical.
    Eigen::VectorXd StartingProbabilitiesInfectedCritical =
        Eigen::VectorXd::Zero((Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedCritical>());
    StartingProbabilitiesInfectedCritical[0] = 0.3;
    StartingProbabilitiesInfectedCritical
        [(Eigen::Index)LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.] = 1 - 0.3;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() = StartingProbabilitiesInfectedCritical;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToDead>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToDead().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2., 7.1);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToRecovered().get_default(
            LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2., 7.1);

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
    Eigen::VectorXd pop(LctState::Count);
    pop[LctState::get_first_index<InfectionState::Susceptible>()] =
        initial_populations[(size_t)InfectionState::Susceptible][0];
    for (size_t i = 0; i < LctState::get_num_subcompartments<InfectionState::Exposed>(); i++) {
        pop[LctState::get_first_index<InfectionState::Exposed>() + i] =
            initial_populations[(size_t)InfectionState::Exposed][i];
    }
    for (size_t i = 0; i < LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>(); i++) {
        pop[LctState::get_first_index<InfectionState::InfectedNoSymptoms>() + i] =
            initial_populations[(size_t)InfectionState::InfectedNoSymptoms][i];
    }
    for (size_t i = 0; i < LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>(); i++) {
        pop[LctState::get_first_index<InfectionState::InfectedSymptoms>() + i] =
            initial_populations[(size_t)InfectionState::InfectedSymptoms][i];
    }
    for (size_t i = 0; i < LctState::get_num_subcompartments<InfectionState::InfectedSevere>(); i++) {
        pop[LctState::get_first_index<InfectionState::InfectedSevere>() + i] =
            initial_populations[(size_t)InfectionState::InfectedSevere][i];
    }
    for (size_t i = 0; i < LctState::get_num_subcompartments<InfectionState::InfectedCritical>(); i++) {
        pop[LctState::get_first_index<InfectionState::InfectedCritical>() + i] =
            initial_populations[(size_t)InfectionState::InfectedCritical][i];
    }
    pop[LctState::get_first_index<InfectionState::Recovered>()] =
        initial_populations[(size_t)InfectionState::Recovered][0];
    pop[LctState::get_first_index<InfectionState::Dead>()] = initial_populations[(size_t)InfectionState::Dead][0];

    // Compare the result of get_derivatives() with a hand calculated result.
    Eigen::VectorXd dydt(LctState::Count);
    model.get_derivatives(pop, pop, 0, dydt);

    Eigen::VectorXd compare(LctState::Count);
    compare << -15.3409, -3.4091, 6.25, -17.5 * 0.91, 15 * 0.91, 0 * 0.91, -17.5 * 0.09, 15 * 0.09, 0 * 0.09,
        3.3052 * 0.2, 3.4483 * 0.2, 3.3052 * 0.8, 3.4483 * 0.8, -7.0417 * 0.25, 6.3158 * 0.25, -7.0417 * 0.75,
        6.3158 * 0.75, -2.2906 * 0.3, -2.8169 * 0.3, -2.2906 * 0.7, -2.8169 * 0.7, 12.3899, 1.6901;

    for (size_t i = 0; i < LctState::Count; i++) {
        ASSERT_NEAR(compare[i], dydt[i], 1e-3) << "Condition failed at index: " << i;
    }
}