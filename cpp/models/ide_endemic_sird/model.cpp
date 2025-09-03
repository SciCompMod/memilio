#include "ide_endemic_sird/model.h"
#include "ide_endemic_sird/computed_parameters.h"
#include "ide_endemic_sird/parameters.h"
#include "ide_endemic_sird/infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/logging.h"

#include "memilio/utils/time_series.h"
#include "vector"
#include <Eigen/src/Core/util/Meta.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

namespace mio
{
namespace endisird
{
Model::Model(CompParameters const& compparams)
    : compparameters{std::make_shared<CompParameters>(compparams)}
    , transitions{TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count))}
    , populations{compparameters->m_statesinit}

{

    // Set flows at start time t0.
    // As we assume that all individuals have infectio age 0 at time t0, the flows at t0 are set to 0.
    transitions.add_time_point(
        0, TimeSeries<ScalarType>::Vector::Constant(static_cast<size_t>(InfectionTransition::Count), 0.));
    // Set population size at start timt t0.
    ScalarType init_populationsize =
        std::accumulate(populations[0].begin(), populations[0].end(), 0) - populations[0][(int)InfectionState::Dead];
    m_totalpopulation.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, init_populationsize));
    m_totalpopulation_derivative.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, 0.));

    // Set the force of infection term at time t0.
    m_forceofinfection.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, compparameters->m_FoI_0[0]));

    //Set normalized_populations at start time t0.
    TimeSeries<ScalarType>::Vector vec_normalizedpopulations =
        TimeSeries<ScalarType>::Vector(Eigen::Index(InfectionState::Count) - 1);
    for (int infection_state = 0; infection_state < Eigen::Index(InfectionState::Count) - 1; infection_state++) {
        vec_normalizedpopulations[infection_state] = populations[0][infection_state] / m_totalpopulation[0][0];
    }
    m_normalizedpopulations.add_time_point(0, vec_normalizedpopulations);
}

bool Model::check_constraints() const
{
    if (!(static_cast<size_t>(populations.get_num_elements()) == static_cast<size_t>(InfectionState::Count))) {
        log_error(" A variable given for model construction is not valid. Number of elements in vector of populations "
                  "does not match the required number.");
        return true;
    }

    for (int i = 0; i < static_cast<int>(InfectionState::Count); i++) {
        if (populations[0][i] < 0) {
            log_error("Initialization failed. Initial values for populations are less than zero.");
            return true;
        }
    }
    return compparameters->check_constraints();
}

// ----Functionality for the iterations of a simulation. ----
void Model::compute_susceptibles(ScalarType dt)
{
    Eigen::Index num_time_points = populations.get_num_time_points();
    populations.get_last_value()[static_cast<int>(InfectionState::Susceptible)] =
        (populations[num_time_points - 2][static_cast<int>(InfectionState::Susceptible)] +
         dt * m_totalpopulation[num_time_points - 2][0] * compparameters->parameters.get<NaturalBirthRate>()) /
        (1 + dt * (m_forceofinfection[num_time_points - 2][0] + compparameters->parameters.get<NaturalDeathRate>()));
}

void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow,
                         Eigen::Index idx_CurrentCompartment, Eigen::Index current_time_index, ScalarType dt)
{
    Eigen::Index calc_time_index =
        (Eigen::Index)std::ceil(compparameters->m_transitiondistributions_support_max[idx_InfectionTransitions] / dt);
    ScalarType current_time_age = static_cast<ScalarType>(current_time_index) * dt;

    ScalarType sum = 0;
    //Determine the starting point of the for loop.
    Eigen::Index starting_point = std::max(0, (int)current_time_index - (int)calc_time_index);

    for (Eigen::Index i = starting_point; i < current_time_index; i++) {
        ScalarType state_age_i = static_cast<ScalarType>(i) * dt;
        sum += transitions[i + 1][idx_IncomingFlow] *
               std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i)) *
               compparameters->m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index - i];
    }
    if (current_time_index <= calc_time_index) {
        transitions.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * compparameters->parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum -
            std::exp((-compparameters->parameters.get<NaturalDeathRate>()) * (current_time_age)) *
                populations[0][idx_CurrentCompartment] *
                compparameters->parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] *
                compparameters->m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index];
    }
    else {
        transitions.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * compparameters->parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum;
    }
}

void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow,
                         Eigen::Index idx_CurrentCompartment, ScalarType dt)
{
    Eigen::Index current_time_index = transitions.get_num_time_points() - 1;
    compute_flow(idx_InfectionTransitions, idx_IncomingFlow, idx_CurrentCompartment, current_time_index, dt);
}

void Model::flows_currents_timestep(ScalarType dt)
{
    Eigen::Index current_time_index = populations.get_num_time_points() - 1;
    // Calculate the transition SusceptibleToExposed with force of infection from previous time step and Susceptibles from
    // current time step.
    transitions.get_last_value()[static_cast<int>(InfectionTransition::SusceptibleToInfected)] =
        m_forceofinfection[current_time_index - 1][0] *
        populations.get_last_value()[static_cast<int>(InfectionState::Susceptible)];

    // Calculate the other Transitions with compute_flow.
    // Infected To Dead:
    compute_flow(Eigen::Index(InfectionTransition::InfectedToDead),
                 Eigen::Index(InfectionTransition::SusceptibleToInfected), Eigen::Index(InfectionState::Infected), dt);
    // Infected To Recovered:
    compute_flow(Eigen::Index(InfectionTransition::InfectedToRecovered),
                 Eigen::Index(InfectionTransition::SusceptibleToInfected), Eigen::Index(InfectionState::Infected), dt);
}

void Model::update_compartment_with_sum(InfectionState infectionState,
                                        std::vector<InfectionTransition> const& IncomingFlows,
                                        bool NaturalDeathispossible, bool Transitionispossible, ScalarType dt)
{
    Eigen::Index current_time_index = populations.get_num_time_points() - 1;
    ScalarType current_time_age     = (ScalarType)current_time_index * dt;
    Eigen::Index calc_time_index    = current_time_index;
    if (Transitionispossible) {
        calc_time_index = compparameters->m_transitiondistributions.size() - 1;
    }

    ScalarType sum = 0;

    Eigen::Index starting_point = std::max(0, (int)current_time_index - (int)calc_time_index);
    for (int i = starting_point; i < current_time_index; i++) {
        ScalarType state_age_i = (ScalarType)i * dt;
        ScalarType sum_inflows = 0;
        for (const InfectionTransition& inflow : IncomingFlows) {
            sum_inflows += transitions[i + 1][(int)inflow];
        }
        if (NaturalDeathispossible && Transitionispossible) {
            sum += compparameters->m_transitiondistributions[current_time_index - i] * sum_inflows *
                   std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i));
        }
        else if (NaturalDeathispossible && !Transitionispossible) {
            sum += sum_inflows *
                   std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i));
        }
        // The case !NaturalDeathispossible && Transitionispossible is not possible, as if NaturalDeath is not possible
        // this means you are in the Death compartment and then Transition is also not possible.
        else {
            sum += sum_inflows;
        }
    }
    if (NaturalDeathispossible && Transitionispossible) {
        if (current_time_index <= calc_time_index) {
            populations.get_last_value()[(int)infectionState] =
                dt * sum + compparameters->m_transitiondistributions[current_time_index] *
                               populations[0][(int)infectionState] *
                               std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age));
        }
        else {
            populations.get_last_value()[(int)infectionState] = dt * sum;
        }
    }
    else if (NaturalDeathispossible && !Transitionispossible) {
        populations.get_last_value()[(int)infectionState] =
            dt * sum + populations[0][(int)infectionState] *
                           std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age));
    }
    else {
        populations.get_last_value()[(int)infectionState] = dt * sum + populations[0][(int)infectionState];
    }
}

void Model::update_compartments(ScalarType dt)
{

    // Infected
    update_compartment_with_sum(InfectionState::Infected, {InfectionTransition::SusceptibleToInfected}, true, true, dt);
    // Recovered
    update_compartment_with_sum(InfectionState::Recovered,
                                {
                                    InfectionTransition::InfectedToRecovered,
                                },
                                true, false, dt);

    // Dead
    update_compartment_with_sum(InfectionState::Dead, {InfectionTransition::InfectedToDead}, false, false, dt);
}

void Model::compute_populationsize()
{
    ScalarType sum = 0;
    for (int state = 0; state < Eigen::Index(InfectionState::Count) - 1; state++) {
        sum += populations.get_last_value()[state];
    }
    m_totalpopulation.get_last_value()[0] = sum;
    // Here we comoute the derivative of the total population that is given by
    // BirthRate * N(t) - DeathRate * N(t) - Transition[InfectedCritical>ToDeath](t).
    m_totalpopulation_derivative.get_last_value()[0] =
        (compparameters->parameters.get<NaturalBirthRate>() - compparameters->parameters.get<NaturalDeathRate>()) *
            m_totalpopulation.get_last_value()[0] -
        transitions.get_last_value()[(int)InfectionTransition::InfectedToDead];
}

void Model::compute_normalizedcompartments()
{
    for (int infection_state = 0; infection_state < Eigen::Index(InfectionState::Count) - 1; infection_state++) {
        m_normalizedpopulations.get_last_value()[infection_state] =
            populations.get_last_value()[infection_state] / m_totalpopulation.get_last_value()[0];
    }
}

void Model::compute_forceofinfection(ScalarType dt)
{

    Eigen::Index num_time_points = populations.get_num_time_points();
    ScalarType current_time      = populations.get_last_time();

    // Determine the starting point of the for loop.
    Eigen::Index starting_point = std::max(0, (int)num_time_points - 1 - (int)compparameters->m_infectivity.size());

    ScalarType sum = 0.;
    // Compute the sum in the force of infection term.
    for (Eigen::Index i = starting_point; i < num_time_points - 1; i++) {
        sum += compparameters->m_infectivity[num_time_points - 1 - i] *
               populations[i + 1][(int)InfectionState::Susceptible] * m_forceofinfection[i][0];
    }

    ScalarType Foi_0 = 0;

    // Add inital functions for the force of infection in case they still have an influence.
    if (num_time_points <= (int)compparameters->m_FoI_0.size()) {
        Foi_0 = compparameters->m_FoI_0[num_time_points - 1] / m_totalpopulation[num_time_points - 1][0];
    }
    m_forceofinfection.get_last_value()[0] =
        (dt * sum * compparameters->parameters.get<TransmissionProbabilityOnContact>().eval(current_time) *
         compparameters->parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(0, 0)) /
            m_totalpopulation[num_time_points - 1][0] +
        Foi_0;
}

}; // namespace endisird

} // namespace mio
