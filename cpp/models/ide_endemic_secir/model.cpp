#include "ide_endemic_secir/model.h"
#include "ide_endemic_secir/computed_parameters.h"
#include "ide_endemic_secir/parameters.h"
#include "ide_endemic_secir/infection_state.h"
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
namespace endisecir
{
Model::Model(CompParameters const& compparams)
    : compparameters{std::make_shared<CompParameters>(compparams)}
    , transitions{TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count))}
    , transitions_update{TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count))}
    , populations{compparameters->m_statesinit}
    , populations_update{compparameters->m_statesinit}

{

    // Set flows at start time t0.
    // As we assume that all individuals have infectio age 0 at time t0, the flows at t0 are set to 0.
    transitions.add_time_point(
        0, TimeSeries<ScalarType>::Vector::Constant(static_cast<size_t>(InfectionTransition::Count), 0.));
    transitions_update.add_time_point(
        0, TimeSeries<ScalarType>::Vector::Constant(static_cast<size_t>(InfectionTransition::Count), 0.));

    // Set population size at start timt t0.
    ScalarType init_populationsize =
        std::accumulate(populations[0].begin(), populations[0].end(), 0) - populations[0][(int)InfectionState::Dead];
    m_totalpopulation.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, init_populationsize));
    m_totalpopulationupdate.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, init_populationsize));

    m_forceofinfection.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, 0));
    m_forceofinfectionupdate.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, 0));

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
// ----Functionality for Initialization. ----
void Model::initialization_compute_forceofinfection()
{
    m_forceofinfection[0][0]       = compparameters->m_FoI_0[0] / m_totalpopulation[0][0];
    m_forceofinfectionupdate[0][0] = compparameters->m_FoI_0[0] / m_totalpopulationupdate[0][0];
}
// ----Functionality for the iterations of a simulation. ----
void Model::compute_susceptibles(ScalarType dt)
{
    Eigen::Index num_time_points = populations.get_num_time_points();
    populations.get_last_value()[static_cast<int>(InfectionState::Susceptible)] =
        (populations[num_time_points - 2][static_cast<int>(InfectionState::Susceptible)] +
         dt * m_totalpopulation[num_time_points - 2][0] * compparameters->parameters.get<NaturalBirthRate>()) /
        (1 + dt * (m_forceofinfection[num_time_points - 2][0] + compparameters->parameters.get<NaturalDeathRate>()));

    //Computation of susceptibles using the update formula for the other compartments.
    //Formula for Susceptibles stays the same, but we use m_totalpopulationupdate and m_forceofinfectionupdate
    populations_update.get_last_value()[static_cast<int>(InfectionState::Susceptible)] =
        (populations_update[num_time_points - 2][static_cast<int>(InfectionState::Susceptible)] +
         dt * (m_totalpopulationupdate[num_time_points - 2][0] * compparameters->parameters.get<NaturalBirthRate>())) /
        (1 +
         dt * (m_forceofinfectionupdate[num_time_points - 2][0] + compparameters->parameters.get<NaturalDeathRate>()));
}

void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow,
                         Eigen::Index idx_CurrentCompartment, Eigen::Index current_time_index, ScalarType dt)
{
    Eigen::Index calc_time_index =
        (Eigen::Index)std::ceil(compparameters->m_transitiondistributions_support_max[idx_InfectionTransitions] / dt);
    ScalarType current_time_age = static_cast<ScalarType>(current_time_index) * dt;

    ScalarType sum1 = 0;
    ScalarType sum2 = 0;
    //Determine the starting point of the for loop.
    Eigen::Index starting_point = std::max(0, (int)current_time_index - (int)calc_time_index);

    for (Eigen::Index i = starting_point; i < current_time_index; i++) {
        ScalarType state_age_i = static_cast<ScalarType>(i) * dt;
        sum1 += transitions[i + 1][idx_IncomingFlow] *
                std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i)) *
                compparameters->m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index - i];
        //For the update formula version of the model:
        sum2 += transitions_update[i + 1][idx_IncomingFlow] *
                std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i)) *
                compparameters->m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index - i];
    }
    if (current_time_index <= calc_time_index) {
        transitions.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * compparameters->parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum1 -
            std::exp((-compparameters->parameters.get<NaturalDeathRate>()) * (current_time_age)) *
                populations[0][idx_CurrentCompartment] *
                compparameters->parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] *
                compparameters->m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index];

        //For the update formula version of the model:
        transitions_update.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * compparameters->parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum2 -
            std::exp((-compparameters->parameters.get<NaturalDeathRate>()) * (current_time_age)) *
                populations_update[0][idx_CurrentCompartment] *
                compparameters->parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] *
                compparameters->m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index];
    }
    else {
        transitions.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * compparameters->parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum1;

        //For the update formula version of the model:
        transitions_update.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * compparameters->parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum2;
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
    transitions.get_last_value()[static_cast<int>(InfectionTransition::SusceptibleToExposed)] =
        m_forceofinfection[current_time_index - 1][0] *
        populations.get_last_value()[static_cast<int>(InfectionState::Susceptible)];

    ///For the update formula version of the model:
    transitions_update.get_last_value()[static_cast<int>(InfectionTransition::SusceptibleToExposed)] =
        m_forceofinfectionupdate[current_time_index - 1][0] *
        populations_update.get_last_value()[static_cast<int>(InfectionState::Susceptible)];

    // Calculate the other Transitions with compute_flow.
    // Exposed to InfectedNoSymptoms
    compute_flow(Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                 Eigen::Index(InfectionTransition::SusceptibleToExposed), Eigen::Index(InfectionState::Exposed), dt);
    // InfectedNoSymptoms to InfectedSymptoms
    compute_flow(Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                 Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                 Eigen::Index(InfectionState::InfectedNoSymptoms), dt);
    // InfectedNoSymptoms to Recovered
    compute_flow(Eigen::Index(InfectionTransition::InfectedNoSymptomsToRecovered),
                 Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                 Eigen::Index(InfectionState::InfectedNoSymptoms), dt);
    // InfectedSymptoms to InfectedSevere
    compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                 Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                 Eigen::Index(InfectionState::InfectedSymptoms), dt);
    // InfectedSymptoms to Recovered
    compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToRecovered),
                 Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                 Eigen::Index(InfectionState::InfectedSymptoms), dt);
    // InfectedSevere to InfectedCritical
    compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                 Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                 Eigen::Index(InfectionState::InfectedSevere), dt);
    // InfectedCritical to Recovered
    compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToRecovered),
                 Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                 Eigen::Index(InfectionState::InfectedSevere), dt);
    // InfectedCritical to Dead
    compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToDead),
                 Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                 Eigen::Index(InfectionState::InfectedCritical), dt);
    // InfectedCritical to Recovered
    compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToRecovered),
                 Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                 Eigen::Index(InfectionState::InfectedCritical), dt);
}

void Model::update_compartment_with_sum(InfectionState infectionState,
                                        std::vector<InfectionTransition> const& IncomingFlows,
                                        bool NaturalDeathispossible, bool Transitionispossible, ScalarType dt)
{
    Eigen::Index current_time_index = populations.get_num_time_points() - 1;
    ScalarType current_time_age     = (ScalarType)current_time_index * dt;
    Eigen::Index calc_time_index    = current_time_index;
    if (Transitionispossible) {
        calc_time_index = compparameters->m_transitiondistributions[(int)infectionState].size() - 1;
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
            sum += compparameters->m_transitiondistributions[(int)infectionState][current_time_index - i] *
                   sum_inflows *
                   std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i));
        }
        else if (NaturalDeathispossible && !Transitionispossible) {
            sum += sum_inflows *
                   std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i));
        }
        else {
            sum += sum_inflows;
        }
    }
    if (NaturalDeathispossible && Transitionispossible) {
        if (current_time_index <= calc_time_index) {
            populations.get_last_value()[(int)infectionState] =
                dt * sum + compparameters->m_transitiondistributions[(int)infectionState][current_time_index] *
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
void Model::update_compartment_from_flow(InfectionState infectionState,
                                         std::vector<InfectionTransition> const& IncomingFlows,
                                         std::vector<InfectionTransition> const& OutgoingFlows,
                                         bool NaturalDeathispossible, ScalarType dt)
{
    Eigen::Index num_time_points = populations.get_num_time_points();

    ScalarType updated_compartment = populations_update[num_time_points - 2][static_cast<int>(infectionState)];

    for (const InfectionTransition& inflow : IncomingFlows) {
        updated_compartment += dt * transitions_update.get_last_value()[(int)inflow];
    }
    for (const InfectionTransition& outflow : OutgoingFlows) {
        updated_compartment -= dt * transitions_update.get_last_value()[(int)outflow];
    }

    if (NaturalDeathispossible) {
        updated_compartment = updated_compartment / (1 + dt * compparameters->parameters.get<NaturalDeathRate>());
    }
    populations_update.get_last_value()[(int)infectionState] = updated_compartment;
}

void Model::update_compartments(ScalarType dt)
{

    // Exposed
    update_compartment_from_flow(InfectionState::Exposed, {InfectionTransition::SusceptibleToExposed},
                                 {InfectionTransition::ExposedToInfectedNoSymptoms}, true, dt);
    update_compartment_with_sum(InfectionState::Exposed, {InfectionTransition::SusceptibleToExposed}, true, true, dt);
    // InfectedNoSymptoms
    update_compartment_from_flow(
        InfectionState::InfectedNoSymptoms, {InfectionTransition::ExposedToInfectedNoSymptoms},
        {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms, InfectionTransition::InfectedNoSymptomsToRecovered},
        true, dt);
    update_compartment_with_sum(InfectionState::InfectedNoSymptoms, {InfectionTransition::ExposedToInfectedNoSymptoms},
                                true, true, dt);

    // InfectedSymptoms
    update_compartment_from_flow(
        InfectionState::InfectedSymptoms, {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms},
        {InfectionTransition::InfectedSymptomsToInfectedSevere, InfectionTransition::InfectedSymptomsToRecovered}, true,
        dt);
    update_compartment_with_sum(InfectionState::InfectedSymptoms,
                                {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms}, true, true, dt);

    // InfectedSevere
    update_compartment_from_flow(
        InfectionState::InfectedSevere, {InfectionTransition::InfectedSymptomsToInfectedSevere},
        {InfectionTransition::InfectedSevereToInfectedCritical, InfectionTransition::InfectedSevereToRecovered}, true,
        dt);
    update_compartment_with_sum(InfectionState::InfectedSevere, {InfectionTransition::InfectedSymptomsToInfectedSevere},
                                true, true, dt);

    // InfectedCritical
    update_compartment_from_flow(
        InfectionState::InfectedCritical, {InfectionTransition::InfectedSevereToInfectedCritical},
        {InfectionTransition::InfectedCriticalToDead, InfectionTransition::InfectedCriticalToRecovered}, true, dt);
    update_compartment_with_sum(InfectionState::InfectedCritical,
                                {InfectionTransition::InfectedSevereToInfectedCritical}, true, true, dt);
    // Recovered
    update_compartment_from_flow(InfectionState::Recovered,
                                 {
                                     InfectionTransition::InfectedNoSymptomsToRecovered,
                                     InfectionTransition::InfectedSymptomsToRecovered,
                                     InfectionTransition::InfectedSevereToRecovered,
                                     InfectionTransition::InfectedCriticalToRecovered,
                                 },
                                 std::vector<InfectionTransition>(), true, dt);
    update_compartment_with_sum(InfectionState::Recovered,
                                {
                                    InfectionTransition::InfectedNoSymptomsToRecovered,
                                    InfectionTransition::InfectedSymptomsToRecovered,
                                    InfectionTransition::InfectedSevereToRecovered,
                                    InfectionTransition::InfectedCriticalToRecovered,
                                },
                                true, false, dt);

    // Dead
    update_compartment_from_flow(InfectionState::Dead, {InfectionTransition::InfectedCriticalToDead},
                                 std::vector<InfectionTransition>(), false, dt);
    update_compartment_with_sum(InfectionState::Dead, {InfectionTransition::InfectedCriticalToDead}, false, false, dt);
}

void Model::compute_populationsize()
{
    ScalarType sum1 = 0;
    ScalarType sum2 = 0;
    for (int state = 0; state < Eigen::Index(InfectionState::Count) - 1; state++) {
        sum1 += populations.get_last_value()[state];
        sum2 += populations_update.get_last_value()[state];
    }
    m_totalpopulation.get_last_value()[0]       = sum1;
    m_totalpopulationupdate.get_last_value()[0] = sum2;
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

    ScalarType sum        = 0.;
    ScalarType sum_update = 0.;
    // Compute the sum in the force of infection term.
    for (Eigen::Index i = starting_point; i < num_time_points - 1; i++) {
        sum += compparameters->m_infectivity[num_time_points - 1 - i] *
               populations[i + 1][(int)InfectionState::Susceptible] * m_forceofinfection[i][0];
        sum_update += compparameters->m_infectivity[num_time_points - 1 - i] *
                      populations_update[i + 1][(int)InfectionState::Susceptible] * m_forceofinfectionupdate[i][0];
    }

    ScalarType Foi_0        = 0;
    ScalarType Foi_0_update = 0;
    ScalarType f            = 0;
    ScalarType f_update     = 0;

    // Add inital functions for the force of infection in case they still have an influence.
    if (num_time_points <= (int)compparameters->m_FoI_0.size()) {
        Foi_0        = compparameters->m_FoI_0[num_time_points - 1] / m_totalpopulation[num_time_points - 1][0];
        Foi_0_update = compparameters->m_FoI_0[num_time_points - 1] / m_totalpopulationupdate[num_time_points - 1][0];
    }
    if (num_time_points <= (int)compparameters->m_InitFoI.size()) {
        f        = compparameters->m_InitFoI[num_time_points - 1] / m_totalpopulation[num_time_points - 1][0];
        f_update = compparameters->m_InitFoI[num_time_points - 1] / m_totalpopulationupdate[num_time_points - 1][0];
    }
    m_forceofinfection.get_last_value()[0] =
        (dt * sum * compparameters->parameters.get<TransmissionProbabilityOnContact>().eval(current_time) *
         compparameters->parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(0, 0)) /
            m_totalpopulation[num_time_points - 1][0] +
        Foi_0 + f;
    m_forceofinfectionupdate.get_last_value()[0] =
        (dt * sum_update *
         (compparameters->parameters.get<TransmissionProbabilityOnContact>().eval(current_time) *
          compparameters->parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(0, 0))) /
            m_totalpopulationupdate[num_time_points - 1][0] +
        Foi_0_update + f_update;
}

}; // namespace endisecir

} // namespace mio
