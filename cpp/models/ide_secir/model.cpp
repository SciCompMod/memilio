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
#include "ide_secir/parameters.h"
#include "infection_state.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"

#include "vector"
#include <algorithm>

namespace mio
{
namespace isecir
{

Model::Model(TimeSeries<ScalarType>&& init, std::vector<ScalarType> N_init, std::vector<ScalarType> deaths,
             int num_agegroups, std::vector<ScalarType> total_confirmed_cases)
    : parameters{Parameters(AgeGroup(num_agegroups))}
    , m_transitions{std::move(init)}
    , m_populations{TimeSeries<ScalarType>(Eigen::Index(InfectionState::Count) * num_agegroups)}
    , m_total_confirmed_cases{total_confirmed_cases}
    , m_N{N_init}
    , m_num_agegroups{num_agegroups}

{
    if (m_transitions.get_num_time_points() > 0) {
        // Add first time point in m_populations according to last time point in m_transitions which is where we start the simulation.
        m_populations.add_time_point<Eigen::VectorXd>(
            m_transitions.get_last_time(),
            TimeSeries<ScalarType>::Vector::Constant((int)InfectionState::Count * m_num_agegroups, 0.));
    }
    else {
        // Initialize m_populations with zero as the first point of time if no data is provided for the transitions.
        // This can happen for example in the case of initialization with real data.
        m_populations.add_time_point<Eigen::VectorXd>(
            0, TimeSeries<ScalarType>::Vector::Constant((int)InfectionState::Count, 0));
    }

    // Set deaths at simulation start time t0.
    for (int group = 0; group < m_num_agegroups; group++) {
        int index                             = get_state_flat_index(Eigen::Index(InfectionState::Dead), group);
        m_populations[Eigen::Index(0)][index] = deaths[static_cast<Eigen::Index>((size_t)group)];
    }
}

// ---- Functionality to calculate the sizes of the compartments for time t0. ----
void Model::compute_compartment_from_flows(ScalarType dt, Eigen::Index idx_InfectionState, AgeGroup group,
                                           Eigen::Index idx_IncomingFlow, int idx_TransitionDistribution1,
                                           int idx_TransitionDistribution2)
{

    ScalarType sum       = 0;
    ScalarType calc_time = 0;
    // Determine relevant calculation area and corresponding index.
    if ((1 - parameters.get<TransitionProbabilities>()[group][idx_TransitionDistribution1]) > 0) {
        calc_time = std::max(
            parameters.get<TransitionDistributions>()[group][idx_TransitionDistribution1].get_support_max(dt, m_tol),
            parameters.get<TransitionDistributions>()[group][idx_TransitionDistribution2].get_support_max(dt, m_tol));
    }
    else {
        calc_time =
            parameters.get<TransitionDistributions>()[group][idx_TransitionDistribution1].get_support_max(dt, m_tol);
    }

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / dt) - 1;

    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {

        ScalarType state_age = (num_time_points - 1 - i) * dt;
        int index            = get_state_flat_index(idx_IncomingFlow, static_cast<Eigen::Index>((size_t)group));

        sum += (parameters.get<TransitionProbabilities>()[group][idx_TransitionDistribution1] *
                    parameters.get<TransitionDistributions>()[group][idx_TransitionDistribution1].eval(state_age) +
                (1 - parameters.get<TransitionProbabilities>()[group][idx_TransitionDistribution1]) *
                    parameters.get<TransitionDistributions>()[group][idx_TransitionDistribution2].eval(state_age)) *
               m_transitions[i + 1][index];
    }

    m_populations.get_last_value()[idx_InfectionState] = sum;
}

void Model::initial_compute_compartments_infection(ScalarType dt)
{
    // We need to compute the initial compartments for every AgeGroup.
    for (auto group = AgeGroup(0); group < AgeGroup(m_num_agegroups); ++group) {
        // Exposed
        compute_compartment_from_flows(dt, Eigen::Index(InfectionState::Exposed), group,
                                       Eigen::Index(InfectionTransition::SusceptibleToExposed),
                                       (int)InfectionTransition::ExposedToInfectedNoSymptoms);
        // InfectedNoSymptoms
        compute_compartment_from_flows(dt, Eigen::Index(InfectionState::InfectedNoSymptoms), group,
                                       Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                                       (int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
                                       (int)InfectionTransition::InfectedNoSymptomsToRecovered);
        // InfectedSymptoms
        compute_compartment_from_flows(dt, Eigen::Index(InfectionState::InfectedSymptoms), group,
                                       Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                                       (int)InfectionTransition::InfectedSymptomsToInfectedSevere,
                                       (int)InfectionTransition::InfectedSymptomsToRecovered);
        // InfectedSevere
        compute_compartment_from_flows(dt, Eigen::Index(InfectionState::InfectedSevere), group,
                                       Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                                       (int)InfectionTransition::InfectedSevereToInfectedCritical,
                                       (int)InfectionTransition::InfectedSevereToRecovered);
        // InfectedCritical
        compute_compartment_from_flows(dt, Eigen::Index(InfectionState::InfectedCritical), group,
                                       Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                                       (int)InfectionTransition::InfectedCriticalToDead,
                                       (int)InfectionTransition::InfectedCriticalToRecovered);
    }
}

void Model::initial_compute_compartments(ScalarType dt)
{
    // The initialization method only affects the Susceptible and Recovered compartments.
    // It is possible to calculate the sizes of the other compartments in advance because only the initial values of the flows are used.
    initial_compute_compartments_infection(dt);

    //Check if there are Susceptibles or Recovered given.
    bool Susceptibles_given = true;
    for (int group = 0; group < m_num_agegroups; ++group) {
        int Si = get_state_flat_index(Eigen::Index(InfectionState::Susceptible), group);
        if (m_populations[Eigen::Index(0)][Si] < 1e-12) {
            Susceptibles_given = false;
            break;
        }
    }
    bool Recovered_given = true;
    for (int group = 0; group < m_num_agegroups; ++group) {
        int Ri = get_state_flat_index(Eigen::Index(InfectionState::Recovered), group);
        if (m_populations[Eigen::Index(0)][Ri] < 1e-12) {
            Recovered_given = false;
            break;
        }
    }

    // We check which Initialization Method we want to use.
    if (!m_total_confirmed_cases.empty() &&
        std::all_of(m_total_confirmed_cases.begin(), m_total_confirmed_cases.end(), [](ScalarType x) {
            return x > 1e-12;
        })) {
        m_initialization_method = 1;
        for (int group = 0; group < m_num_agegroups; ++group) {
            int Si = get_state_flat_index(Eigen::Index(InfectionState::Susceptible), group);
            int Ei = get_state_flat_index(Eigen::Index(InfectionState::Exposed), group);
            int Ci = get_state_flat_index(Eigen::Index(InfectionState::InfectedNoSymptoms), group);
            int Ii = get_state_flat_index(Eigen::Index(InfectionState::InfectedSymptoms), group);
            int Hi = get_state_flat_index(Eigen::Index(InfectionState::InfectedSevere), group);
            int Ui = get_state_flat_index(Eigen::Index(InfectionState::InfectedCritical), group);
            int Ri = get_state_flat_index(Eigen::Index(InfectionState::Recovered), group);
            int Di = get_state_flat_index(Eigen::Index(InfectionState::Dead), group);
            // The scheme of the ODE model for initialization is applied here.
            m_populations[Eigen::Index(0)][Ri] = m_total_confirmed_cases[group] - m_populations[Eigen::Index(0)][Ii] -
                                                 m_populations[Eigen::Index(0)][Hi] -
                                                 m_populations[Eigen::Index(0)][Ui] -
                                                 m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Dead)];

            m_populations[Eigen::Index(0)][Si] =
                m_N[group] - m_populations[Eigen::Index(0)][Ei] - m_populations[Eigen::Index(0)][Ci] -
                m_populations[Eigen::Index(0)][Ii] - m_populations[Eigen::Index(0)][Hi] -
                m_populations[Eigen::Index(0)][Ui] - m_populations[Eigen::Index(0)][Ri] -
                m_populations[Eigen::Index(0)][Di];
        }
    }

    else if (Susceptibles_given) {
        // Take initialized value for Susceptibles if value can't be calculated via the standard formula.
        m_initialization_method = 2;

        // R; need an initial value for R, therefore do not calculate via compute_recovered()
        for (int group = 0; group < m_num_agegroups; ++group) {
            int Si = get_state_flat_index(Eigen::Index(InfectionState::Susceptible), group);
            int Ei = get_state_flat_index(Eigen::Index(InfectionState::Exposed), group);
            int Ci = get_state_flat_index(Eigen::Index(InfectionState::InfectedNoSymptoms), group);
            int Ii = get_state_flat_index(Eigen::Index(InfectionState::InfectedSymptoms), group);
            int Hi = get_state_flat_index(Eigen::Index(InfectionState::InfectedSevere), group);
            int Ui = get_state_flat_index(Eigen::Index(InfectionState::InfectedCritical), group);
            int Ri = get_state_flat_index(Eigen::Index(InfectionState::Recovered), group);
            int Di = get_state_flat_index(Eigen::Index(InfectionState::Dead), group);
            m_populations[Eigen::Index(0)][Ri] =
                m_N[group] - m_populations[Eigen::Index(0)][Si] - m_populations[Eigen::Index(0)][Ei] -
                m_populations[Eigen::Index(0)][Ci] - m_populations[Eigen::Index(0)][Ii] -
                m_populations[Eigen::Index(0)][Hi] - m_populations[Eigen::Index(0)][Ui] -
                m_populations[Eigen::Index(0)][Di];
        }
    }
    else if (Recovered_given) {
        // If value for Recovered is initialized and standard method is not applicable (i.e., no value for total infections
        // or Susceptibles given directly), calculate Susceptibles via other compartments.
        m_initialization_method = 3;
        for (int group = 0; group < m_num_agegroups; ++group) {
            int Si = get_state_flat_index(Eigen::Index(InfectionState::Susceptible), group);
            int Ei = get_state_flat_index(Eigen::Index(InfectionState::Exposed), group);
            int Ci = get_state_flat_index(Eigen::Index(InfectionState::InfectedNoSymptoms), group);
            int Ii = get_state_flat_index(Eigen::Index(InfectionState::InfectedSymptoms), group);
            int Hi = get_state_flat_index(Eigen::Index(InfectionState::InfectedSevere), group);
            int Ui = get_state_flat_index(Eigen::Index(InfectionState::InfectedCritical), group);
            int Ri = get_state_flat_index(Eigen::Index(InfectionState::Recovered), group);
            int Di = get_state_flat_index(Eigen::Index(InfectionState::Dead), group);
            m_populations[Eigen::Index(0)][Si] =
                m_N[group] - m_populations[Eigen::Index(0)][Ei] - m_populations[Eigen::Index(0)][Ci] -
                m_populations[Eigen::Index(0)][Ii] - m_populations[Eigen::Index(0)][Hi] -
                m_populations[Eigen::Index(0)][Ui] - m_populations[Eigen::Index(0)][Ri] -
                m_populations[Eigen::Index(0)][Di];
        }
    }
    else {
        // Compute Susceptibles at t0 and m_forceofinfection at time t0-dt as initial values for discretization scheme.
        // Use m_forceofinfection at t0-dt to be consistent with further calculations of S (see compute_susceptibles()),
        // where also the value of m_forceofinfection for the previous timestep is used.
        compute_forceofinfection(dt, true);
        if (std::all_of(m_forceofinfection.begin(), m_forceofinfection.end(), [](ScalarType x) {
                return x > 1e-12;
            })) {
            m_initialization_method = 4;
            for (int group = 0; group < m_num_agegroups; ++group) {
                int Si = get_state_flat_index(Eigen::Index(InfectionState::Susceptible), group);
                int Ei = get_state_flat_index(Eigen::Index(InfectionState::Exposed), group);
                int Ci = get_state_flat_index(Eigen::Index(InfectionState::InfectedNoSymptoms), group);
                int Ii = get_state_flat_index(Eigen::Index(InfectionState::InfectedSymptoms), group);
                int Hi = get_state_flat_index(Eigen::Index(InfectionState::InfectedSevere), group);
                int Ui = get_state_flat_index(Eigen::Index(InfectionState::InfectedCritical), group);
                int Ri = get_state_flat_index(Eigen::Index(InfectionState::Recovered), group);
                int Di = get_state_flat_index(Eigen::Index(InfectionState::Dead), group);

                int SEi = get_transition_flat_index(Eigen::Index(InfectionTransition::SusceptibleToExposed), group);
                /* Attention: With an inappropriate combination of parameters and initial conditions, it can happen that S 
                becomes greater than N when this formula is used. In this case, at least one compartment is initialized 
                with a number less than zero, so that a log_error is thrown.
                However, this initialization method is consistent with the numerical solver of the model equations,
                so it may sometimes make sense to use this method. */
                m_populations[Eigen::Index(0)][Si] =
                    m_transitions.get_last_value()[SEi] / (dt * m_forceofinfection[group]);

                // Recovered; calculated as the difference between the total population and the sum of the other compartment sizes.
                m_populations[Eigen::Index(0)][Ri] =
                    m_N[group] - m_populations[Eigen::Index(0)][Si] - m_populations[Eigen::Index(0)][Ei] -
                    m_populations[Eigen::Index(0)][Ci] - m_populations[Eigen::Index(0)][Ii] -
                    m_populations[Eigen::Index(0)][Hi] - m_populations[Eigen::Index(0)][Ui] -
                    m_populations[Eigen::Index(0)][Di];
            }
        }
        else {
            m_initialization_method = -1;
            log_error("Error occured while initializing compartments: Force of infection is evaluated to 0 and neither "
                      "Susceptibles nor Recovered or total confirmed cases for time 0 were set. One of them should be "
                      "larger 0.");
        }
    }
    // Verify that the initialization produces an appropriate result.
    // Another check would be if the sum of the compartments is equal to N, but in all initialization methods one of the compartments is initialized via N - the others.
    // This also means that if a compartment is greater than N, we will always have one or more compartments less than zero.
    // Check if all compartments are non negative.
    for (int group = 0; group < m_num_agegroups; ++group) {
        for (int i = 0; i < (int)InfectionState::Count; i++) {
            int idx = get_state_flat_index(i, group);
            if (m_populations[0][idx] < 0) {
                m_initialization_method = -2;
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
    // we need to compute to compute the Susceptibles in every AgeGroup.
    for (int group = 0; group < m_num_agegroups; ++group) {
        Eigen::Index num_time_points = m_populations.get_num_time_points();
        // Using number of Susceptibles from previous time step and force of infection from previous time step:
        // Compute current number of Susceptibles and store Susceptibles in m_populations.
        int Si = get_state_flat_index(Eigen::Index(InfectionState::Susceptible), group);
        m_populations.get_last_value()[Si] =
            m_populations[num_time_points - 2][Si] / (1 + dt * m_forceofinfection[group]);
    }
}

void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt,
                         Eigen::Index current_time_index, AgeGroup group)
{
    ScalarType sum = 0;
    /* In order to satisfy TransitionDistribution(dt*i) = 0 for all i >= k, k is determined by the maximum support of the distribution.
    Since we are using a backwards difference scheme to compute the derivative, we have that the
    derivative of TransitionDistribution(dt*i) = 0 for all i >= k+1.

    Hence calc_time_index goes until std::ceil(support_max/dt) since for std::ceil(support_max/dt)+1 all terms are already zero. 
    This needs to be adjusted if we are changing the finite difference scheme */

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(
        parameters.get<TransitionDistributions>()[group][idx_InfectionTransitions].get_support_max(dt, m_tol) / dt);

    int Transition_idx = get_transition_flat_index(idx_InfectionTransitions, static_cast<Eigen::Index>((size_t)group));
    for (Eigen::Index i = current_time_index - calc_time_index; i < current_time_index; i++) {
        // (current_time_index - i) * dt is the time the individuals have already spent in this state.
        ScalarType state_age = (current_time_index - i) * dt;

        // Use backward difference scheme to approximate first derivative.
        sum += (parameters.get<TransitionDistributions>()[group][idx_InfectionTransitions].eval(state_age) -
                parameters.get<TransitionDistributions>()[group][idx_InfectionTransitions].eval(state_age - dt)) /
               dt * m_transitions[i + 1][idx_IncomingFlow];
    }

    m_transitions.get_value(current_time_index)[Transition_idx] =
        (-dt) * parameters.get<TransitionProbabilities>()[group][idx_InfectionTransitions] * sum;
}

void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt,
                         AgeGroup group)
{
    Eigen::Index current_time_index = m_transitions.get_num_time_points() - 1;
    compute_flow(idx_InfectionTransitions, idx_IncomingFlow, dt, current_time_index, group);
}

void Model::flows_current_timestep(ScalarType dt)
{
    //We need to compute the flows for every AgeGroup.
    for (auto group = AgeGroup(0); group < AgeGroup(m_num_agegroups); ++group) {
        int SEi = get_transition_flat_index(Eigen::Index(InfectionTransition::SusceptibleToExposed),
                                            static_cast<Eigen::Index>((size_t)group));
        int Si =
            get_state_flat_index(Eigen::Index(InfectionState::Susceptible), static_cast<Eigen::Index>((size_t)group));

        // Calculate flow SusceptibleToExposed with force of infection from previous time step and Susceptibles from current time step.
        m_transitions.get_last_value()[SEi] =
            dt * m_forceofinfection[static_cast<Eigen::Index>((size_t)group)] * m_populations.get_last_value()[Si];

        // Calculate all other flows with compute_flow.
        // Flow from Exposed to InfectedNoSymptoms
        compute_flow(Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                     Eigen::Index(InfectionTransition::SusceptibleToExposed), dt, group);
        // Flow from InfectedNoSymptoms to InfectedSymptoms
        compute_flow(Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                     Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, group);
        // Flow from InfectedNoSymptoms to Recovered
        compute_flow(Eigen::Index(InfectionTransition::InfectedNoSymptomsToRecovered),
                     Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, group);
        // Flow from InfectedSymptoms to InfectedSevere
        compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                     Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt, group);
        // Flow from InfectedSymptoms to Recovered
        compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToRecovered),
                     Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt, group);
        // Flow from InfectedSevere to InfectedCritical
        compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                     Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, group);
        // Flow from InfectedSevere to Recovered
        compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToRecovered),
                     Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, group);
        // Flow from InfectedCritical to Dead
        compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToDead),
                     Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, group);
        // Flow from InfectedCritical to Recovered
        compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToRecovered),
                     Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, group);
    }
}

void Model::update_compartments()
{
    //We need to update the compartments for every AgeGroup.
    for (auto group = AgeGroup(0); group < AgeGroup(m_num_agegroups); ++group) {
        // Exposed
        update_compartment_from_flow(InfectionState::Exposed, {InfectionTransition::SusceptibleToExposed},
                                     {InfectionTransition::ExposedToInfectedNoSymptoms}, group);

        // InfectedNoSymptoms
        update_compartment_from_flow(InfectionState::InfectedNoSymptoms,
                                     {InfectionTransition::ExposedToInfectedNoSymptoms},
                                     {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
                                      InfectionTransition::InfectedNoSymptomsToRecovered},
                                     group);

        // InfectedSymptoms
        update_compartment_from_flow(
            InfectionState::InfectedSymptoms, {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms},
            {InfectionTransition::InfectedSymptomsToInfectedSevere, InfectionTransition::InfectedSymptomsToRecovered},
            group);

        // InfectedSevere
        update_compartment_from_flow(
            InfectionState::InfectedSevere, {InfectionTransition::InfectedSymptomsToInfectedSevere},
            {InfectionTransition::InfectedSevereToInfectedCritical, InfectionTransition::InfectedSevereToRecovered},
            group);

        // InfectedCritical
        update_compartment_from_flow(
            InfectionState::InfectedCritical, {InfectionTransition::InfectedSevereToInfectedCritical},
            {InfectionTransition::InfectedCriticalToDead, InfectionTransition::InfectedCriticalToRecovered}, group);

        // Recovered
        update_compartment_from_flow(InfectionState::Recovered,
                                     {
                                         InfectionTransition::InfectedNoSymptomsToRecovered,
                                         InfectionTransition::InfectedSymptomsToRecovered,
                                         InfectionTransition::InfectedSevereToRecovered,
                                         InfectionTransition::InfectedCriticalToRecovered,
                                     },
                                     std::vector<InfectionTransition>(), group);

        // Dead
        update_compartment_from_flow(InfectionState::Dead, {InfectionTransition::InfectedCriticalToDead},
                                     std::vector<InfectionTransition>(), group);
    }
}

void Model::update_compartment_from_flow(InfectionState infectionState,
                                         std::vector<InfectionTransition> const& IncomingFlows,
                                         std::vector<InfectionTransition> const& OutgoingFlows, AgeGroup group)
{
    int state_idx = get_state_flat_index(Eigen::Index(infectionState), static_cast<Eigen::Index>((size_t)group));

    Eigen::Index num_time_points   = m_populations.get_num_time_points();
    ScalarType updated_compartment = m_populations[num_time_points - 2][state_idx];
    for (const InfectionTransition& inflow : IncomingFlows) {
        int inflow_idx = get_transition_flat_index(Eigen::Index(inflow), static_cast<Eigen::Index>((size_t)group));
        updated_compartment += m_transitions.get_last_value()[inflow_idx];
    }
    for (const InfectionTransition& outflow : OutgoingFlows) {
        int outflow_idx = get_transition_flat_index(Eigen::Index(outflow), static_cast<Eigen::Index>((size_t)group));
        updated_compartment -= m_transitions.get_last_value()[outflow_idx];
    }
    m_populations.get_last_value()[state_idx] = updated_compartment;
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

        int Di  = get_state_flat_index(Eigen::Index(InfectionState::Dead), static_cast<Eigen::Index>((size_t)i));
        int HDi = get_transition_flat_index(Eigen::Index(InfectionTransition::InfectedCriticalToDead),
                                            static_cast<Eigen::Index>((size_t)i));

        if (initialization) {
            // Determine m_forceofinfection at time t0-dt which is the penultimate timepoint in m_transitions.
            num_time_points = m_transitions.get_num_time_points() - 1;
            // Get time of penultimate timepoint in m_transitions.
            current_time = m_transitions.get_time(num_time_points - 1);

            // Determine the number of individuals in Dead compartment at time t0-dt.
            deaths_i = m_populations[Eigen::Index(0)][Di] - m_transitions.get_last_value()[HDi];
        }
        else {
            // Determine m_forceofinfection for current last time in m_transitions.
            num_time_points = m_transitions.get_num_time_points();
            current_time    = m_transitions.get_last_time();
            deaths_i        = m_populations.get_last_value()[Di];
        }
        // We need to sum, over contacts with all Age Groups.
        for (auto j = AgeGroup(0); j < AgeGroup(m_num_agegroups); ++j) {

            ScalarType sum = 0;

            for (Eigen::Index l = num_time_points - 1 - calc_time_index; l < num_time_points - 1; l++) {

                int ECj = get_transition_flat_index(Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                                                    static_cast<Eigen::Index>((size_t)j));
                int CRj = get_transition_flat_index(Eigen::Index(InfectionTransition::InfectedNoSymptomsToRecovered),
                                                    static_cast<Eigen::Index>((size_t)j));

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
                         m_transitions[l + 1][ECj] *
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
                         m_transitions[l + 1][CRj] *
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
    for (auto group = AgeGroup(0); group < AgeGroup(m_num_agegroups); ++group) {
        global_support_max_new = std::max(
            {parameters.get<TransitionDistributions>()[group][(int)InfectionTransition::ExposedToInfectedNoSymptoms]
                 .get_support_max(dt, m_tol),
             parameters
                 .get<TransitionDistributions>()[group][(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[group][(int)InfectionTransition::InfectedNoSymptomsToRecovered]
                 .get_support_max(dt, m_tol),
             parameters
                 .get<TransitionDistributions>()[group][(int)InfectionTransition::InfectedSymptomsToInfectedSevere]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[group][(int)InfectionTransition::InfectedSymptomsToRecovered]
                 .get_support_max(dt, m_tol),
             parameters
                 .get<TransitionDistributions>()[group][(int)InfectionTransition::InfectedSevereToInfectedCritical]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[group][(int)InfectionTransition::InfectedSevereToRecovered]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[group][(int)InfectionTransition::InfectedCriticalToDead]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[group][(int)InfectionTransition::InfectedCriticalToRecovered]
                 .get_support_max(dt, m_tol)});
        if (global_support_max_new > global_support_max) {
            global_support_max = global_support_max_new;
        }
    }
    return global_support_max;
}

} // namespace isecir
} // namespace mio
