/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Anna Wendler, Lena Ploetzke, Martin J. Kuehn
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
#include "Eigen/src/Core/util/Meta.h"
#include "ide_secir/parameters.h"
#include "infection_state.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"

#include "vector"

namespace mio
{
namespace isecir
{

Model::Model(CustomIndexArray<TimeSeries<ScalarType>, AgeGroup>&& init, std::vector<ScalarType> N_init,
             std::vector<ScalarType> deaths, int num_agegroups, std::vector<ScalarType> total_confirmed_cases,
             const ParameterSet& Parameterset_init)
    : parameters{Parameterset_init}
    , m_transitions{std::move(init)}
    , m_populations{CustomIndexArray<TimeSeries<ScalarType>, AgeGroup>(AgeGroup(num_agegroups),
                                                                       Eigen::Index(InfectionState::Count))}
    , m_total_confirmed_cases{total_confirmed_cases}
    , m_N{N_init}
    , m_num_agegroups{num_agegroups}

{
    for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); ++i) {
        if (m_transitions[i].get_num_time_points() > 0) {
            // Add first time point in m_populations according to last time point in m_transitions which is where we start the simulation.
            m_populations[i].add_time_point<Eigen::VectorXd>(
                m_transitions[i].get_last_time(),
                TimeSeries<ScalarType>::Vector::Constant((int)InfectionState::Count, 0.));
        }
        else {
            // Initialize m_populations with zero as the first point of time if no data is provided for the transitions.
            // This can happen for example in the case of initialization with real data.
            m_populations[i].add_time_point<Eigen::VectorXd>(
                0, TimeSeries<ScalarType>::Vector::Constant((int)InfectionState::Count, 0));
        }

        // Set deaths at simulation start time t0.
        m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Dead)] =
            deaths[static_cast<Eigen::Index>((size_t)i)];
    }
}

// ---- Functionality to calculate the sizes of the compartments for time t0. ----
void Model::compute_compartment_from_flows(ScalarType dt, Eigen::Index idx_InfectionState,
                                           Eigen::Index idx_IncomingFlow, int idx_TransitionDistribution1,
                                           int idx_TransitionDistribution2)
{
    for (auto j = AgeGroup(0); j < AgeGroup(m_num_agegroups); ++j) {
        ScalarType sum       = 0;
        ScalarType calc_time = 0;
        // Determine relevant calculation area and corresponding index.
        if ((1 - parameters.get<TransitionProbabilities>()[j][idx_TransitionDistribution1]) > 0) {
            calc_time = std::max(
                parameters.get<TransitionDistributions>()[j][idx_TransitionDistribution1].get_support_max(dt, m_tol),
                parameters.get<TransitionDistributions>()[j][idx_TransitionDistribution2].get_support_max(dt, m_tol));
        }
        else {
            calc_time =
                parameters.get<TransitionDistributions>()[j][idx_TransitionDistribution1].get_support_max(dt, m_tol);
        }

        Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / dt) - 1;

        Eigen::Index num_time_points = m_transitions[j].get_num_time_points();

        for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {

            ScalarType state_age = (num_time_points - 1 - i) * dt;

            sum += (parameters.get<TransitionProbabilities>()[j][idx_TransitionDistribution1] *
                        parameters.get<TransitionDistributions>()[j][idx_TransitionDistribution1].eval(state_age) +
                    (1 - parameters.get<TransitionProbabilities>()[j][idx_TransitionDistribution1]) *
                        parameters.get<TransitionDistributions>()[j][idx_TransitionDistribution2].eval(state_age)) *
                   m_transitions[j][i + 1][idx_IncomingFlow];
        }

        m_populations[j].get_last_value()[idx_InfectionState] = sum;
    }
}

void Model::initial_compute_compartments_infection(ScalarType dt)
{
    // Exposed
    compute_compartment_from_flows(dt, Eigen::Index(InfectionState::Exposed),
                                   Eigen::Index(InfectionTransition::SusceptibleToExposed),
                                   (int)InfectionTransition::ExposedToInfectedNoSymptoms);
    // InfectedNoSymptoms
    compute_compartment_from_flows(dt, Eigen::Index(InfectionState::InfectedNoSymptoms),
                                   Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                                   (int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
                                   (int)InfectionTransition::InfectedNoSymptomsToRecovered);
    // InfectedSymptoms
    compute_compartment_from_flows(dt, Eigen::Index(InfectionState::InfectedSymptoms),
                                   Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                                   (int)InfectionTransition::InfectedSymptomsToInfectedSevere,
                                   (int)InfectionTransition::InfectedSymptomsToRecovered);
    // InfectedSevere
    compute_compartment_from_flows(dt, Eigen::Index(InfectionState::InfectedSevere),
                                   Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                                   (int)InfectionTransition::InfectedSevereToInfectedCritical,
                                   (int)InfectionTransition::InfectedSevereToRecovered);
    // InfectedCritical
    compute_compartment_from_flows(dt, Eigen::Index(InfectionState::InfectedCritical),
                                   Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                                   (int)InfectionTransition::InfectedCriticalToDead,
                                   (int)InfectionTransition::InfectedCriticalToRecovered);
}

void Model::initial_compute_compartments(ScalarType dt)
{
    // The initialization method only affects the Susceptible and Recovered compartments.
    // It is possible to calculate the sizes of the other compartments in advance because only the initial values of the flows are used.
    initial_compute_compartments_infection(dt);
    for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); ++i) {
        if (m_total_confirmed_cases[static_cast<Eigen::Index>((size_t)i)] > 1e-12) {
            m_initialization_method[static_cast<Eigen::Index>((size_t)i)] = 1;

            // The scheme of the ODE model for initialization is applied here.
            m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] =
                m_total_confirmed_cases[static_cast<Eigen::Index>((size_t)i)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSymptoms)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSevere)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedCritical)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Dead)];

            m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] =
                m_N[static_cast<Eigen::Index>((size_t)i)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Exposed)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedNoSymptoms)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSymptoms)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSevere)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedCritical)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Dead)];
        }
        else if (m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] > 1e-12) {
            // Take initialized value for Susceptibles if value can't be calculated via the standard formula.
            m_initialization_method[static_cast<Eigen::Index>((size_t)i)] = 2;

            // R; need an initial value for R, therefore do not calculate via compute_recovered()
            m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] =
                m_N[static_cast<Eigen::Index>((size_t)i)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Exposed)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedNoSymptoms)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSymptoms)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSevere)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedCritical)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Dead)];
        }
        else if (m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] > 1e-12) {
            // If value for Recovered is initialized and standard method is not applicable (i.e., no value for total infections
            // or Susceptibles given directly), calculate Susceptibles via other compartments.
            m_initialization_method[static_cast<Eigen::Index>((size_t)i)] = 3;

            m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] =
                m_N[static_cast<Eigen::Index>((size_t)i)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Exposed)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedNoSymptoms)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSymptoms)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSevere)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedCritical)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] -
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Dead)];
        }
        else {
            // Compute Susceptibles at t0 and m_forceofinfection at time t0-dt as initial values for discretization scheme.
            // Use m_forceofinfection at t0-dt to be consistent with further calculations of S (see compute_susceptibles()),
            // where also the value of m_forceofinfection for the previous timestep is used.
            compute_forceofinfection(dt, true);
            if (m_forceofinfection[static_cast<Eigen::Index>((size_t)i)] > 1e-12) {
                m_initialization_method[static_cast<Eigen::Index>((size_t)i)] = 4;

                /* Attention: With an inappropriate combination of parameters and initial conditions, it can happen that S 
                becomes greater than N when this formula is used. In this case, at least one compartment is initialized 
                with a number less than zero, so that a log_error is thrown.
                However, this initialization method is consistent with the numerical solver of the model equations,
                so it may sometimes make sense to use this method. */
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] =
                    m_transitions[i].get_last_value()[Eigen::Index(InfectionTransition::SusceptibleToExposed)] /
                    (dt * m_forceofinfection[static_cast<Eigen::Index>((size_t)i)]);

                // Recovered; calculated as the difference between the total population and the sum of the other compartment sizes.
                m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] =
                    m_N[static_cast<Eigen::Index>((size_t)i)] -
                    m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] -
                    m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Exposed)] -
                    m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedNoSymptoms)] -
                    m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSymptoms)] -
                    m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSevere)] -
                    m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::InfectedCritical)] -
                    m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Dead)];
            }
            else {
                m_initialization_method[static_cast<Eigen::Index>((size_t)i)] = -1;
                log_error(
                    "Error occured while initializing compartments: Force of infection is evaluated to 0 and neither "
                    "Susceptibles nor Recovered or total confirmed cases for time 0 were set. One of them should be "
                    "larger 0.");
            }
        }
        // Verify that the initialization produces an appropriate result.
        // Another check would be if the sum of the compartments is equal to N, but in all initialization methods one of the compartments is initialized via N - the others.
        // This also means that if a compartment is greater than N, we will always have one or more compartments less than zero.
        // Check if all compartments are non negative.
        for (int j = 0; j < (int)InfectionState::Count; j++) {
            if (m_populations[i][0][j] < 0) {
                m_initialization_method[static_cast<Eigen::Index>((size_t)i)] = -2;
                log_error(
                    "Initialization failed. One or more initial values for populations are less than zero. It may "
                    "help to use a different method for initialization.");
            }
        }
    }

    // Compute m_forceofinfection at time t0 needed for further simulation.
    compute_forceofinfection(dt);
}

// ---- Functionality for the iterations of a simulation. ----
void Model::compute_susceptibles(ScalarType dt)
{
    for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); ++i) {
        Eigen::Index num_time_points = m_populations[i].get_num_time_points();
        // Using number of Susceptibles from previous time step and force of infection from previous time step:
        // Compute current number of Susceptibles and store Susceptibles in m_populations.
        m_populations[i].get_last_value()[Eigen::Index(InfectionState::Susceptible)] =
            m_populations[i][num_time_points - 2][Eigen::Index(InfectionState::Susceptible)] /
            (1 + dt * m_forceofinfection[static_cast<Eigen::Index>((size_t)i)]);
    }
}

void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt,
                         Eigen::Index current_time_index, AgeGroup agegroup)
{
    ScalarType sum = 0;
    /* In order to satisfy TransitionDistribution(dt*i) = 0 for all i >= k, k is determined by the maximum support of the distribution.
    Since we are using a backwards difference scheme to compute the derivative, we have that the
    derivative of TransitionDistribution(dt*i) = 0 for all i >= k+1.

    Hence calc_time_index goes until std::ceil(support_max/dt) since for std::ceil(support_max/dt)+1 all terms are already zero. 
    This needs to be adjusted if we are changing the finite difference scheme */

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(
        parameters.get<TransitionDistributions>()[agegroup][idx_InfectionTransitions].get_support_max(dt, m_tol) / dt);

    for (Eigen::Index i = current_time_index - calc_time_index; i < current_time_index; i++) {
        // (current_time_index - i) * dt is the time the individuals have already spent in this state.
        ScalarType state_age = (current_time_index - i) * dt;

        // Use backward difference scheme to approximate first derivative.
        sum += (parameters.get<TransitionDistributions>()[agegroup][idx_InfectionTransitions].eval(state_age) -
                parameters.get<TransitionDistributions>()[agegroup][idx_InfectionTransitions].eval(state_age - dt)) /
               dt * m_transitions[agegroup][i + 1][idx_IncomingFlow];
    }

    m_transitions[agegroup].get_value(current_time_index)[Eigen::Index(idx_InfectionTransitions)] =
        (-dt) * parameters.get<TransitionProbabilities>()[agegroup][idx_InfectionTransitions] * sum;
}

void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt,
                         AgeGroup i)
{
    Eigen::Index current_time_index = m_transitions[i].get_num_time_points() - 1;
    compute_flow(idx_InfectionTransitions, idx_IncomingFlow, dt, current_time_index, i);
}

void Model::flows_current_timestep(ScalarType dt)
{
    for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); ++i) {
        // Calculate flow SusceptibleToExposed with force of infection from previous time step and Susceptibles from current time step.
        m_transitions[i].get_last_value()[Eigen::Index(InfectionTransition::SusceptibleToExposed)] =
            dt * m_forceofinfection[static_cast<Eigen::Index>((size_t)i)] *
            m_populations[i].get_last_value()[Eigen::Index(InfectionState::Susceptible)];

        // Calculate all other flows with compute_flow.
        // Flow from Exposed to InfectedNoSymptoms
        compute_flow(Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                     Eigen::Index(InfectionTransition::SusceptibleToExposed), dt, i);
        // Flow from InfectedNoSymptoms to InfectedSymptoms
        compute_flow(Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                     Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, i);
        // Flow from InfectedNoSymptoms to Recovered
        compute_flow(Eigen::Index(InfectionTransition::InfectedNoSymptomsToRecovered),
                     Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, i);
        // Flow from InfectedSymptoms to InfectedSevere
        compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                     Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt, i);
        // Flow from InfectedSymptoms to Recovered
        compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToRecovered),
                     Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt, i);
        // Flow from InfectedSevere to InfectedCritical
        compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                     Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, i);
        // Flow from InfectedSevere to Recovered
        compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToRecovered),
                     Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, i);
        // Flow from InfectedCritical to Dead
        compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToDead),
                     Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, i);
        // Flow from InfectedCritical to Recovered
        compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToRecovered),
                     Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, i);
    }
}

void Model::update_compartments()
{
    for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); ++i) {
        // Exposed
        update_compartment_from_flow(InfectionState::Exposed, {InfectionTransition::SusceptibleToExposed},
                                     {InfectionTransition::ExposedToInfectedNoSymptoms}, i);

        // InfectedNoSymptoms
        update_compartment_from_flow(InfectionState::InfectedNoSymptoms,
                                     {InfectionTransition::ExposedToInfectedNoSymptoms},
                                     {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
                                      InfectionTransition::InfectedNoSymptomsToRecovered},
                                     i);

        // InfectedSymptoms
        update_compartment_from_flow(
            InfectionState::InfectedSymptoms, {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms},
            {InfectionTransition::InfectedSymptomsToInfectedSevere, InfectionTransition::InfectedSymptomsToRecovered},
            i);

        // InfectedSevere
        update_compartment_from_flow(
            InfectionState::InfectedSevere, {InfectionTransition::InfectedSymptomsToInfectedSevere},
            {InfectionTransition::InfectedSevereToInfectedCritical, InfectionTransition::InfectedSevereToRecovered}, i);

        // InfectedCritical
        update_compartment_from_flow(
            InfectionState::InfectedCritical, {InfectionTransition::InfectedSevereToInfectedCritical},
            {InfectionTransition::InfectedCriticalToDead, InfectionTransition::InfectedCriticalToRecovered}, i);

        // Recovered
        update_compartment_from_flow(InfectionState::Recovered,
                                     {
                                         InfectionTransition::InfectedNoSymptomsToRecovered,
                                         InfectionTransition::InfectedSymptomsToRecovered,
                                         InfectionTransition::InfectedSevereToRecovered,
                                         InfectionTransition::InfectedCriticalToRecovered,
                                     },
                                     std::vector<InfectionTransition>(), i);

        // Dead
        update_compartment_from_flow(InfectionState::Dead, {InfectionTransition::InfectedCriticalToDead},
                                     std::vector<InfectionTransition>(), i);
    }
}

void Model::update_compartment_from_flow(InfectionState infectionState,
                                         std::vector<InfectionTransition> const& IncomingFlows,
                                         std::vector<InfectionTransition> const& OutgoingFlows, AgeGroup agegroup)
{
    Eigen::Index num_time_points   = m_populations[agegroup].get_num_time_points();
    ScalarType updated_compartment = m_populations[agegroup][num_time_points - 2][Eigen::Index(infectionState)];
    for (const InfectionTransition& inflow : IncomingFlows) {
        updated_compartment += m_transitions[agegroup].get_last_value()[Eigen::Index(inflow)];
    }
    for (const InfectionTransition& outflow : OutgoingFlows) {
        updated_compartment -= m_transitions[agegroup].get_last_value()[Eigen::Index(outflow)];
    }
    m_populations[agegroup].get_last_value()[Eigen::Index(infectionState)] = updated_compartment;
}

void Model::compute_forceofinfection(ScalarType dt, bool initialization)
{

    m_forceofinfection = std::vector<ScalarType>(m_num_agegroups, 0);
    // We compute the force of infection for every AgeGroup.
    for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); ++i) {

        // Determine the relevant calculation area = union of the supports of the relevant transition distributions.
        ScalarType calc_time = std::max(
            {parameters
                 .get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedNoSymptomsToRecovered]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedSymptomsToInfectedSevere]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedSymptomsToRecovered]
                 .get_support_max(dt, m_tol)});

        // Corresponding index.
        // Need calc_time_index timesteps in sum,
        // subtract 1 because in the last summand all TransitionDistributions evaluate to 0 (by definition of support_max).
        Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / dt) - 1;

        Eigen::Index num_time_points;
        ScalarType current_time;
        ScalarType deaths_i;

        if (initialization) {
            // Determine m_forceofinfection at time t0-dt which is the penultimate timepoint in m_transitions.
            num_time_points = m_transitions[i].get_num_time_points() - 1;
            // Get time of penultimate timepoint in m_transitions.
            current_time = m_transitions[i].get_time(num_time_points - 1);

            // Determine the number of individuals in Dead compartment at time t0-dt.
            deaths_i = m_populations[i][Eigen::Index(0)][Eigen::Index(InfectionState::Dead)] -
                       m_transitions[i].get_last_value()[Eigen::Index(InfectionTransition::InfectedCriticalToDead)];
        }
        else {
            // Determine m_forceofinfection for current last time in m_transitions.
            num_time_points = m_transitions[i].get_num_time_points();
            current_time    = m_transitions[i].get_last_time();
            deaths_i        = m_populations[i].get_last_value()[Eigen::Index(InfectionState::Dead)];
        }
        // We need to sum, over contacts with all Age Groups.
        for (auto j = AgeGroup(0); j < AgeGroup(m_num_agegroups); ++j) {

            ScalarType sum = 0;

            for (Eigen::Index l = num_time_points - 1 - calc_time_index; l < num_time_points - 1; l++) {

                ScalarType state_age = (num_time_points - 1 - l) * dt;
                ScalarType season_val =
                    1 + parameters.get<Seasonality>() *
                            sin(3.141592653589793 *
                                (std::fmod((parameters.get<StartDay>() + current_time), 365.0) / 182.5 + 0.5));
                sum +=
                    season_val * parameters.get<TransmissionProbabilityOnContact>()[i].eval(state_age) *
                    parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(
                        static_cast<Eigen::Index>((size_t)i), static_cast<Eigen::Index>((size_t)j)) *
                    ((parameters.get<TransitionProbabilities>()[j][(
                          int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
                          parameters
                              .get<TransitionDistributions>()[j][(
                                  int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
                              .eval(state_age) +
                      parameters.get<TransitionProbabilities>()[j][(
                          int)InfectionTransition::InfectedNoSymptomsToRecovered] *
                          parameters
                              .get<TransitionDistributions>()[j]
                                                             [(int)InfectionTransition::InfectedNoSymptomsToRecovered]
                              .eval(state_age)) *
                         m_transitions[j][l + 1][Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms)] *
                         parameters.get<RelativeTransmissionNoSymptoms>()[j].eval(state_age) +
                     (parameters.get<TransitionProbabilities>()[j][(
                          int)InfectionTransition::InfectedSymptomsToInfectedSevere] *
                          parameters
                              .get<TransitionDistributions>()[j][(
                                  int)InfectionTransition::InfectedSymptomsToInfectedSevere]
                              .eval(state_age) +
                      parameters.get<TransitionProbabilities>()[j]
                                                               [(int)InfectionTransition::InfectedSymptomsToRecovered] *
                          parameters
                              .get<TransitionDistributions>()[j][(int)InfectionTransition::InfectedSymptomsToRecovered]
                              .eval(state_age)) *
                         m_transitions[j][l + 1]
                                      [Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] *
                         parameters.get<RiskOfInfectionFromSymptomatic>()[j].eval(state_age));
            }
            m_forceofinfection[static_cast<Eigen::Index>((size_t)i)] +=
                1 / (m_N[static_cast<Eigen::Index>((size_t)i)] - deaths_i) * sum;
        }
    }
}

ScalarType Model::get_global_support_max(ScalarType dt) const
{
    ScalarType global_support_max     = 0.;
    ScalarType global_support_max_new = 0.;
    for (auto i = AgeGroup(0); i < AgeGroup(m_num_agegroups); ++i) {
        global_support_max_new = std::max(
            {parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::ExposedToInfectedNoSymptoms]
                 .get_support_max(dt, m_tol),
             parameters
                 .get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedNoSymptomsToRecovered]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedSymptomsToInfectedSevere]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedSymptomsToRecovered]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedSevereToInfectedCritical]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedSevereToRecovered]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedCriticalToDead]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[i][(int)InfectionTransition::InfectedCriticalToRecovered]
                 .get_support_max(dt, m_tol)});
        if (global_support_max_new > global_support_max) {
            global_support_max = global_support_max_new;
        }
    }
    return global_support_max;
}

} // namespace isecir
} // namespace mio
