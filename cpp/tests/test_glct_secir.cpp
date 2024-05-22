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
#include "glct_secir/simulation.h"
#include "glct_secir/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"

#include <gtest/gtest.h>

// Test if the function eval_right_hand_side() is working using a hand calculated result.
TEST(TestGLCTSecir, testEvalRightHandSide)
{
    // Define model.
    // Chose more than one subcompartment for all compartments except S, R, D so that the function is correct for all selections.
    using Model    = mio::glsecir::Model<2, 3, 2, 2, 2>;
    using LctState = Model::LctState;

    // Define initial population distribution in infection states, one entry per subcompartment.
    Eigen::VectorXd init(LctState::Count);
    init[LctState::get_first_index<LctState::InfectionState::Susceptible>()]            = 750;
    init[LctState::get_first_index<LctState::InfectionState::Exposed>()]                = 30;
    init[LctState::get_first_index<LctState::InfectionState::Exposed>() + 1]            = 20;
    init[LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>()]     = 20;
    init[LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>() + 1] = 10;
    init[LctState::get_first_index<LctState::InfectionState::InfectedNoSymptoms>() + 2] = 10;
    init[LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>()]       = 30;
    init[LctState::get_first_index<LctState::InfectionState::InfectedSymptoms>() + 1]   = 20;
    init[LctState::get_first_index<LctState::InfectionState::InfectedSevere>()]         = 40;
    init[LctState::get_first_index<LctState::InfectionState::InfectedSevere>() + 1]     = 10;
    init[LctState::get_first_index<LctState::InfectionState::InfectedCritical>()]       = 10;
    init[LctState::get_first_index<LctState::InfectionState::InfectedCritical>() + 1]   = 20;
    init[LctState::get_first_index<LctState::InfectionState::Recovered>()]              = 20;
    init[LctState::get_first_index<LctState::InfectionState::Dead>()]                   = 10;

    Model model(std::move(init));

    // Set parameters.
    // Exposed.
    model.parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::Exposed>());
    model.parameters.get<mio::glsecir::TransitionMatrixExposed>() = mio::glsecir::TransitionMatrixExposed().get_default(
        LctState::get_num_subcompartments<LctState::InfectionState::Exposed>(), 3.2);
    // InfectedNoSymptoms.
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        mio::glsecir::StartingProbabilitiesInfectedNoSymptoms().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedNoSymptoms>());
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptoms>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptoms().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedNoSymptoms>(), 2.);
    // InfectedSymptoms.
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() =
        mio::glsecir::StartingProbabilitiesInfectedSymptoms().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedSymptoms>());
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptoms>() =
        mio::glsecir::TransitionMatrixInfectedSymptoms().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedSymptoms>(), 5.8);
    // InfectedSevere.
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() =
        mio::glsecir::StartingProbabilitiesInfectedSevere().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedSevere>());
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevere>() =
        mio::glsecir::TransitionMatrixInfectedSevere().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedSevere>(), 9.5);
    // InfectedCritical.
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() =
        mio::glsecir::StartingProbabilitiesInfectedCritical().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedCritical>());
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCritical>() =
        mio::glsecir::TransitionMatrixInfectedCritical().get_default(
            LctState::get_num_subcompartments<LctState::InfectionState::InfectedCritical>(), 7.1);

    model.parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::glsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));

    model.parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model.parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    model.parameters.get<mio::glsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    model.parameters.get<mio::glsecir::Seasonality>()                    = 0.;
    model.parameters.get<mio::glsecir::StartDay>()                       = 0;
    model.parameters.get<mio::glsecir::SeverePerInfectedSymptoms>()      = 0.2;
    model.parameters.get<mio::glsecir::CriticalPerSevere>()              = 0.25;
    model.parameters.get<mio::glsecir::DeathsPerCritical>()              = 0.3;

    // Compare the result of eval_right_hand_side() with a hand calculated result.
    int num_subcompartments = LctState::Count;
    Eigen::VectorXd dydt(num_subcompartments);
    model.eval_right_hand_side(model.get_initial_values(), 0, dydt);

    Eigen::VectorXd compare(num_subcompartments);
    compare << -15.3409, -3.4091, 6.25, -17.5, 15, 0, 3.3052, 3.4483, -7.0417, 6.3158, -2.2906, -2.8169, 12.3899,
        1.6901;

    for (int i = 0; i < num_subcompartments; i++) {
        ASSERT_NEAR(compare[i], dydt[i], 1e-3) << "Condition failed at index: " << i;
    }
}