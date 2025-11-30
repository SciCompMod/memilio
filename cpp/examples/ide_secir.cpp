/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler, Lena Ploetzke, Hannah Tritzschak
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
    // This is a simple example to demonstrate how use the IDE-SECIR model.

    using Vec = Eigen::VectorX<ScalarType>;

    // Define simulation parameters.
    ScalarType t0   = 0.;
    ScalarType tmax = 5.;
    ScalarType dt   = 0.01; // The step size will stay constant throughout the simulation.

    // Define number of age groups.
    size_t num_agegroups = 2;

    // Define initial values for the total population and number of deaths per age group.
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population_init =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 1000.);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths_init =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 6.);

    // Create TimeSeries with num_transitions * num_agegroups elements where initial transitions needed for simulation
    // will be stored. We require values for the transitions for a sufficient number of time points before the start of
    // the simulation to initialize our model.
    size_t num_transitions = (size_t)mio::isecir::InfectionTransition::Count;
    mio::TimeSeries<ScalarType> transitions_init(num_transitions * num_agegroups);

    // Define vector of transitions that will be added as values to the time points of the TimeSeries transitions_init.
    Vec vec_init(num_transitions * num_agegroups);
    for (size_t group = 0; group < num_agegroups; ++group) {
        vec_init[group * num_transitions + (size_t)mio::isecir::InfectionTransition::SusceptibleToExposed] = 25.0;
        vec_init[group * num_transitions + (size_t)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] =
            15.0;
        vec_init[group * num_transitions +
                 (size_t)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
        vec_init[group * num_transitions + (size_t)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
            4.0;
        vec_init[group * num_transitions + (size_t)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
            1.0;
        vec_init[group * num_transitions + (size_t)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered] = 4.0;
        vec_init[group * num_transitions + (size_t)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] =
            1.0;
        vec_init[group * num_transitions + (size_t)mio::isecir::InfectionTransition::InfectedSevereToRecovered]   = 1.0;
        vec_init[group * num_transitions + (size_t)mio::isecir::InfectionTransition::InfectedCriticalToDead]      = 1.0;
        vec_init[group * num_transitions + (size_t)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] = 1.0;
    }
    // Multiply vec_init with dt so that within a time interval of length 1, always the above number of
    // individuals are transitioning from one compartment to another, irrespective of the chosen time step size.
    vec_init = vec_init * dt;

    // In this example, we will set the TransitionDistributions below. For these distributions, setting the initial time
    // point of the TimeSeries transitions_init at time -10 will give us a sufficient number of time points before t0=0.
    // For more information on this, we refer to the documentation of TransitionDistributions in
    // models/ide_secir/parameters.h.
    transitions_init.add_time_point(-10, vec_init);
    // Add further time points with distance dt until time t0.
    while (transitions_init.get_last_time() < t0 - dt / 2) {
        transitions_init.add_time_point(transitions_init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(transitions_init), total_population_init, deaths_init, num_agegroups);

    // Uncomment one of the code blocks below to use a different method to initialize the model, based on a
    // given number of either Susceptibles or Recovered instead of using the TimeSeries transitions_init from above.

    // Initialization method with given Susceptibles.
    // size_t num_infstates = (size_t)mio::isecir::InfectionState::Count;
    // for (size_t group = 0; group < num_agegroups; ++group) {
    //     model.populations.get_last_value()[group * num_infstates + (size_t)mio::isecir::InfectionState::Susceptible] =
    //         900;
    // }

    // Initialization method with given Recovered.
    // size_t num_infstates = (size_t)mio::isecir::InfectionState::Count;
    // for (size_t group = 0; group < num_agegroups; ++group) {
    //     model.populations.get_last_value()[group * num_infstates + (size_t)mio::isecir::InfectionState::Recovered] = 10;
    // }

    // Set working parameters.

    // TransitionDistributions
    // In the following, we explicitly set the TransitionDistributions for the first age group. If the model contains
    // more age groups, the default distributions are used for these age groups.
    mio::SmootherCosine<ScalarType> smoothcos1(3.0);
    mio::StateAgeFunctionWrapper<ScalarType> delaydistribution1(smoothcos1);
    std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_delaydistrib1(num_transitions, delaydistribution1);
    // TransitionDistribution is not used for SusceptibleToExposed. Therefore, the parameter can be set to any value.
    vec_delaydistrib1[(size_t)mio::isecir::InfectionTransition::SusceptibleToExposed].set_distribution_parameter(-1.);
    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(0)] = vec_delaydistrib1;

    // TransitionProbabilities
    std::vector<ScalarType> vec_prob(num_transitions, 0.5);
    // The following probabilities must be 1, as there is no other way to go.
    vec_prob[(size_t)mio::isecir::InfectionTransition::SusceptibleToExposed]        = 1;
    vec_prob[(size_t)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] = 1;
    for (mio::AgeGroup group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        model.parameters.get<mio::isecir::TransitionProbabilities>()[group] = vec_prob;
    }

    // Contact patterns
    mio::ContactMatrixGroup<ScalarType> contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, num_agegroups);
    contact_matrix[0] =
        mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(num_agegroups, num_agegroups, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // Furhter epidemiological parameters
    mio::ExponentialSurvivalFunction<ScalarType> exponential(0.5);
    mio::StateAgeFunctionWrapper<ScalarType> prob(exponential);
    for (mio::AgeGroup group = mio::AgeGroup(0); group < mio::AgeGroup(num_agegroups); ++group) {
        model.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[group] = prob;
        model.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[group]   = prob;
        model.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[group]   = prob;
    }
    model.parameters.set<mio::isecir::Seasonality>(0.1);
    model.parameters.set<mio::isecir::StartDay>(
        40); // Start the simulation on the 40th day of a year (i.e. in February).

    // Check if all model constraints regarding initial values and parameters are satisfied before simulating.
    model.check_constraints(dt);

    // Carry out simulation.
    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);

    // Interpolate results to days.
    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result(), dt / 2.);

    // Print results. Note that the column labels are suitable for a simulation with two age groups and may need to be
    // adapted when the number of age groups is changed.
    // interpolated_results.print_table(
    //     {"S1", "E1", "C1", "I1", "H1", "U1", "R1", "D1 ", "S2", "E2", "C2", "I2", "H2", "U2", "R2", "D2 "}, 16, 8);
    // Uncomment this line to print the transitions.
    // sim.get_transitions().print_table({"S->E 1", "E->C 1", "C->I 1", "C->R 1", "I->H 1", "I->R 1", "H->U 1",
    //                                    "H->R 1", "U->D 1", "U->R 1", "S->E 2", "E->C 2", "C->I 2", "C->R 2",
    //                                    "I->H 2", "I->R 2", "H->U 2", "H->R 2", "U->D 2", "U->R 2"},
    //                                   16, 8);
}
