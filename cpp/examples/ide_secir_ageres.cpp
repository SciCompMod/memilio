/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Hannah Tritzschak, Anna Wendler
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
#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/data/analyze_result.h"

int main()
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    size_t num_agegroups = 2;

    ScalarType tmax = 5;
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> N =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 10000.);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 6.);
    ScalarType dt = 1.;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions * num_agegroups elements where transitions needed for simulation will be
    // stored.
    mio::TimeSeries<ScalarType> init(num_transitions * num_agegroups);

    // Add time points for initialization of transitions.
    Vec vec_init(num_transitions * num_agegroups);
    // Values for the Infectiontransitions are the same for all AgeGroups.
    for (size_t group = 0; group < num_agegroups; ++group) {
        vec_init[group * num_transitions + (int)mio::isecir::InfectionTransition::SusceptibleToExposed]          = 25.0;
        vec_init[group * num_transitions + (int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]   = 15.0;
        vec_init[group * num_transitions +
                 (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]                    = 8.0;
        vec_init[group * num_transitions + (int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] = 4.0;
        vec_init[group * num_transitions + (int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
            1.0;
        vec_init[group * num_transitions + (int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered] = 4.0;
        vec_init[group * num_transitions + (int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] =
            1.0;
        vec_init[group * num_transitions + (int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]   = 1.0;
        vec_init[group * num_transitions + (int)mio::isecir::InfectionTransition::InfectedCriticalToDead]      = 1.0;
        vec_init[group * num_transitions + (int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] = 1.0;
    }

    // Add initial time point to time series.
    init.add_time_point(-10, vec_init);
    // Add further time points until time 0.
    while (init.get_last_time() < -dt / 2) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, deaths, num_agegroups);

    // Uncomment these lines to use a different method to initialize the model using the TimeSeries init.
    // Initialization method with Susceptibles.
    // model.populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Susceptible] = 1000;
    // model.populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Count +
    //                                      (Eigen::Index)mio::isecir::InfectionState::Susceptible] = 1000;
    // Initialization method with Recovered.
    // model.populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered] = 0;
    // model.populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Count +
    //                                      (Eigen::Index)mio::isecir::InfectionState::Recovered] = 0;

    // Set working parameters.
    // First AgeGroup for Transition Distributions.
    mio::SmootherCosine smoothcos1(2.0);
    mio::StateAgeFunctionWrapper delaydistribution1(smoothcos1);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib1(num_transitions, delaydistribution1);
    // TransitionDistribution is not used for SusceptibleToExposed. Therefore, the parameter can be set to any value.
    vec_delaydistrib1[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_distribution_parameter(-1.);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(0)] = vec_delaydistrib1;

    //Second AgeGroup for Transition Distributions.
    mio::SmootherCosine smoothcos2(3.0);
    mio::StateAgeFunctionWrapper delaydistribution2(smoothcos2);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib2(num_transitions, delaydistribution2);
    // TransitionDistribution is not used for SusceptibleToExposed. Therefore, the parameter can be set to any value.
    vec_delaydistrib2[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_distribution_parameter(-1.);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(1)] = vec_delaydistrib2;

    std::vector<ScalarType> vec_prob(num_transitions, 0.5);
    // The following probabilities must be 1, as there is no other way to go.
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    for (mio::AgeGroup group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        model.parameters.get<mio::isecir::TransitionProbabilities>()[group] = vec_prob;
    }

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, static_cast<Eigen::Index>(num_agegroups));
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_agegroups, num_agegroups, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ExponentialSurvivalFunction exponential(0.5);
    mio::StateAgeFunctionWrapper prob(exponential);
    for (mio::AgeGroup group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        model.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[group] = prob;
        model.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[group]   = prob;
        model.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[group]   = prob;
    }
    model.parameters.set<mio::isecir::Seasonality>(0.1);
    // Start the simulation on the 40th day of a year (i.e. in February).
    model.parameters.set<mio::isecir::StartDay>(40);

    model.check_constraints(dt);

    // Carry out simulation.
    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result(), dt / 2.);

    interpolated_results.print_table(
        {"S1", "E1", "C1", "I1", "H1", "U1", "R1", "D1 ", "S2", "E2", "C2", "I2", "H2", "U2", "R2", "D2 "}, 16, 8);
    // Uncomment this line to print the transitions.
    // sim.get_transitions().print_table({"S->E 1", "E->C 1", "C->I 1", "C->R 1", "I->H 1", "I->R 1", "H->U 1",
    //                                    "H->R 1", "U->D 1", "U->R 1", "S->E 2", "E->C 2", "C->I 2", "C->R 2",
    //                                    "I->H 2", "I->R 2", "H->U 2", "H->R 2", "U->D 2", "U->R 2"},
    //                                   16, 8);
}
