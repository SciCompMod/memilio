#include "ide_endemic_secir/model.h"
#include "ide_endemic_secir/parameters.h"
#include "ide_endemic_secir/infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/logging.h"

#include "memilio/utils/time_series.h"
#include "ode_secir/parameters.h"
#include "vector"
#include <Eigen/src/Core/util/Meta.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

namespace mio
{
namespace endisecir
{
Model::Model(TimeSeries<ScalarType>&& states_init)
    : parameters{Parameters()}
    , transitions{TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count))}
    , transitions_update{TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count))}
    , populations{std::move(states_init)}
    , populations_update{TimeSeries<ScalarType>(Eigen::Index(InfectionState::Count))}
{
    populations_update.add_time_point(0, populations[0]);
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
    // m_totalpopulationincludingD.add_time_point(
    //     0,
    //     TimeSeries<ScalarType>::Vector::Constant(1, init_populationsize + populations[0][(int)InfectionState::Dead]));
    // m_totalpopulationupdateincludingD.add_time_point(
    //     0,
    //     TimeSeries<ScalarType>::Vector::Constant(1, init_populationsize + populations[0][(int)InfectionState::Dead]));
    // Set the force of infection at start time t0.
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
    return parameters.check_constraints();
}

// ----Functionality for the iterations of a simulation. ----
void Model::compute_susceptibles(ScalarType dt)
{
    Eigen::Index num_time_points = populations.get_num_time_points();
    populations.get_last_value()[static_cast<int>(InfectionState::Susceptible)] =
        (populations[num_time_points - 2][static_cast<int>(InfectionState::Susceptible)] *
             (1 - dt * parameters.get<NaturalDeathRate>()) +
         dt * m_totalpopulation[num_time_points - 2][0] * parameters.get<NaturalBirthRate>()) /
        (1 + dt * m_forceofinfection[num_time_points - 2][0]);

    //Computation of susceptibles using the update formula for the other compartments.
    //Formula for Susceptibles stays the same, but we use m_totalpopulationupdate and m_forceofinfectionupdate
    populations_update.get_last_value()[static_cast<int>(InfectionState::Susceptible)] =
        (populations_update[num_time_points - 2][static_cast<int>(InfectionState::Susceptible)] *
             (1 - dt * parameters.get<NaturalDeathRate>()) +
         dt * m_totalpopulationupdate[num_time_points - 2][0] * parameters.get<NaturalBirthRate>()) /
        (1 + dt * m_forceofinfectionupdate[num_time_points - 2][0]);
}

void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow,
                         Eigen::Index idx_CurrentCompartment, Eigen::Index current_time_index, ScalarType dt)
{
    Eigen::Index calc_time_index =
        (Eigen::Index)std::ceil(m_transitiondistributions_support_max[idx_InfectionTransitions] / dt);
    ScalarType current_time_age = static_cast<ScalarType>(current_time_index) * dt;

    ScalarType sum1 = 0;
    ScalarType sum2 = 0;
    //Determine the starting point of the for loop.
    Eigen::Index starting_point = std::max(0, (int)current_time_index - (int)calc_time_index);

    for (Eigen::Index i = starting_point; i < current_time_index; i++) {
        ScalarType state_age_i = static_cast<ScalarType>(i) * dt;
        sum1 += transitions[i + 1][idx_IncomingFlow] *
                std::exp(-parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i)) *
                m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index - i];
        //For the update formula version of the model:
        sum2 += transitions_update[i + 1][idx_IncomingFlow] *
                std::exp(-parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i)) *
                m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index - i];
    }
    if (current_time_index <= calc_time_index) {
        transitions.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum1 -
            std::exp((-parameters.get<NaturalDeathRate>()) * (current_time_age)) *
                populations[0][idx_CurrentCompartment] *
                parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] *
                m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index];

        //For the update formula version of the model:
        transitions_update.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum2 -
            std::exp((-parameters.get<NaturalDeathRate>()) * (current_time_age)) *
                populations_update[0][idx_CurrentCompartment] *
                parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] *
                m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index];
    }
    else {
        transitions.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum1;

        //For the update formula version of the model:
        transitions_update.get_value(current_time_index)[idx_InfectionTransitions] =
            (-dt) * parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum2;
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
    ScalarType sum                  = 0;

    for (int i = 0; i < current_time_index; i++) {
        ScalarType state_age_i = (ScalarType)i * dt;
        ScalarType sum_inflows = 0;
        for (const InfectionTransition& inflow : IncomingFlows) {
            sum_inflows += transitions[i + 1][(int)inflow];
        }
        if (NaturalDeathispossible && Transitionispossible) {
            sum += m_transitiondistributions[(int)infectionState][current_time_index - i] * sum_inflows *
                   std::exp(-parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i));
        }
        else if (NaturalDeathispossible && !Transitionispossible) {
            sum += sum_inflows * std::exp(-parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i));
        }
        else {
            sum += sum_inflows;
        }
    }
    if (NaturalDeathispossible && Transitionispossible) {
        populations.get_last_value()[(int)infectionState] =
            dt * sum + m_transitiondistributions[(int)infectionState][current_time_index] *
                           populations[0][(int)infectionState] *
                           std::exp(-parameters.get<NaturalDeathRate>() * (current_time_age));
    }
    else if (NaturalDeathispossible && !Transitionispossible) {
        populations.get_last_value()[(int)infectionState] =
            dt * sum +
            populations[0][(int)infectionState] * std::exp(-parameters.get<NaturalDeathRate>() * (current_time_age));
    }
    else {
        populations.get_last_value()[(int)infectionState] =
            dt * sum +
            m_transitiondistributions[(int)infectionState][current_time_index] * populations[0][(int)infectionState];
    }
}
void Model::update_compartment_from_flow(InfectionState infectionState,
                                         std::vector<InfectionTransition> const& IncomingFlows,
                                         std::vector<InfectionTransition> const& OutgoingFlows,
                                         bool NaturalDeathispossible, ScalarType dt)
{
    Eigen::Index num_time_points = populations.get_num_time_points();

    ScalarType updated_compartment = 0;
    if (NaturalDeathispossible) {
        updated_compartment = populations_update[num_time_points - 2][static_cast<int>(infectionState)] *
                              (1 - dt * parameters.get<NaturalDeathRate>());
    }
    else {
        updated_compartment = populations_update[num_time_points - 2][static_cast<int>(infectionState)];
    }
    for (const InfectionTransition& inflow : IncomingFlows) {
        updated_compartment += dt * transitions_update.get_last_value()[(int)inflow];
    }
    for (const InfectionTransition& outflow : OutgoingFlows) {
        updated_compartment -= dt * transitions_update.get_last_value()[(int)outflow];
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
    m_totalpopulation.get_last_value()[0] =
        std::accumulate(populations.get_last_value().begin(), populations.get_last_value().end(), 0) -
        populations.get_last_value()[(int)InfectionState::Dead];
    m_totalpopulationupdate.get_last_value()[0] =
        std::accumulate(populations_update.get_last_value().begin(), populations_update.get_last_value().end(), 0) -
        populations_update.get_last_value()[(int)InfectionState::Dead];
    // m_totalpopulationincludingD.get_last_value()[0] =
    //     std::accumulate(populations.get_last_value().begin(), populations.get_last_value().end(), 0);
    // m_totalpopulationupdateincludingD.get_last_value()[0] =
    //     std::accumulate(populations_update.get_last_value().begin(), populations_update.get_last_value().end(), 0);
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

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(m_calctime / dt) - 1;

    ScalarType foi_0 = 0.;
    ScalarType f     = 0.;

    //Otherwise the initial values do not have an influence anymore.
    if (num_time_points - 1 <= calc_time_index) {

        //Computation of phi_0:
        foi_0 =
            m_meaninfectivity[num_time_points - 1] *
            (populations[0][static_cast<int>(InfectionState::InfectedNoSymptoms)] *
                 m_transitiondistributions[static_cast<int>(InfectionState::InfectedNoSymptoms)][num_time_points - 1] +
             populations[0][static_cast<int>(InfectionState::InfectedSymptoms)] *
                 m_transitiondistributions[static_cast<int>(InfectionState::InfectedSymptoms)][num_time_points - 1]);

        //Computation of f:
        ScalarType f_1 = 0.;
        ScalarType f_2 = 0.;
        ScalarType f_3 = 0.;

        for (int i = 1; i < num_time_points; i++) {
            ScalarType state_age_i = static_cast<ScalarType>(i - 1) * dt;
            //Compute sum for f_1:
            f_1 +=
                m_transitiondistributions_derivative[static_cast<int>(InfectionTransition::ExposedToInfectedNoSymptoms)]
                                                    [num_time_points - 1 - i] *
                parameters.get<RelativeTransmissionNoSymptoms>().eval(state_age_i) *
                m_transitiondistributions[static_cast<int>(InfectionState::InfectedNoSymptoms)][i - 1];

            //Compute sum for f_2:
            f_2 += m_transitiondistributions_derivative[static_cast<int>(
                       InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)][num_time_points - 1 - i] *
                   parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age_i) *
                   m_transitiondistributions[static_cast<int>(InfectionState::InfectedSymptoms)][i - 1];

            //Compute sum for f_3:
            ScalarType sum_f = 0.;
            for (int j = 1; j <= i; j++) {
                ScalarType state_age_j = static_cast<ScalarType>(j - 1) * dt;
                sum_f += parameters.get<TransitionProbabilities>()[static_cast<int>(
                             InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] *
                         m_transitiondistributions_derivative[static_cast<int>(
                             InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)][num_time_points - 1 - j] *
                         parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age_j) *
                         m_transitiondistributions[static_cast<int>(InfectionState::InfectedSymptoms)][j - 1];
            }
            f_3 +=
                m_transitiondistributions_derivative[static_cast<int>(InfectionTransition::ExposedToInfectedNoSymptoms)]
                                                    [num_time_points - 1 - i] *
                sum_f;
        }
        f = (dt * dt * populations[0][static_cast<int>(InfectionState::Exposed)] *
                 std::exp(-parameters.get<NaturalBirthRate>() * current_time) * f_3 -
             dt * populations[0][static_cast<int>(InfectionState::InfectedNoSymptoms)] *
                 std::exp(-parameters.get<NaturalBirthRate>() * current_time) * f_2 -
             dt * populations[0][static_cast<int>(InfectionState::Exposed)] *
                 std::exp(-parameters.get<NaturalBirthRate>() * current_time) * f_1);
    }

    // Computation of the rest

    // Determine the starting point of the for loop.
    Eigen::Index starting_point = std::max(0, (int)num_time_points - 1 - (int)calc_time_index);

    ScalarType sum = 0.;
    for (Eigen::Index i = starting_point; i < num_time_points - 1; i++) {
        sum += m_meaninfectivity[num_time_points - 1 - i] * populations[i + 1][(int)InfectionState::Susceptible] *
               m_forceofinfection[i][0];
    }

    m_forceofinfection.get_last_value()[0] =
        foi_0 / m_totalpopulation[num_time_points - 1][0] + f / m_totalpopulation[num_time_points - 1][0] +
        (dt * sum) *
            (parameters.get<TransmissionProbabilityOnContact>().eval(current_time) *
             parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(0, 0)) /
            m_totalpopulation[num_time_points - 1][0];

    //For the update formula version of the model:
    sum = 0.;
    for (Eigen::Index i = starting_point; i < num_time_points - 1; i++) {
        sum += m_meaninfectivity[num_time_points - 1 - i] *
               populations_update[i + 1][(int)InfectionState::Susceptible] * m_forceofinfectionupdate[i][0];
    }

    m_forceofinfectionupdate.get_last_value()[0] =
        foi_0 / m_totalpopulationupdate[num_time_points - 1][0] + f / m_totalpopulationupdate[num_time_points - 1][0] +
        (dt * sum) *
            (parameters.get<TransmissionProbabilityOnContact>().eval(current_time) *
             parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(0, 0)) /
            m_totalpopulationupdate[num_time_points - 1][0];
}

// ---- Functionality to set vectors with necessary information regarding TransitionDistributions. ----
void Model::set_transitiondistributions_support_max(ScalarType dt)
{
    // The transition Susceptible to Exposed is not needed in the computations.
    for (int transition = 1; transition < (int)InfectionTransition::Count; transition++) {
        m_transitiondistributions_support_max[transition] =
            parameters.get<TransitionDistributions>()[transition].get_support_max(dt, m_tol);
    }
}

void Model::set_calctime()
{
    // Transitions needed in the force of infection term.
    std::vector<std::vector<int>> relevant_transitions = {
        {(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
         (int)InfectionTransition::InfectedNoSymptomsToRecovered},
        {(int)InfectionTransition::InfectedSymptomsToInfectedSevere,
         (int)InfectionTransition::InfectedSymptomsToRecovered}};

    m_calctime = std::max({m_transitiondistributions_support_max[relevant_transitions[0][0]],
                           m_transitiondistributions_support_max[relevant_transitions[0][1]],
                           m_transitiondistributions_support_max[relevant_transitions[1][0]],
                           m_transitiondistributions_support_max[relevant_transitions[1][1]]});
}

void Model::set_transitiondistributions(ScalarType dt)
{
    //Exposed state:
    Eigen::Index support_max_index = (Eigen::Index)std::ceil(
        m_transitiondistributions_support_max[(int)InfectionTransition::ExposedToInfectedNoSymptoms] / dt);

    std::vector<ScalarType> vec_contribution_to_foi_1(support_max_index + 1, 0.);
    for (Eigen::Index i = 0; i <= support_max_index; i++) {
        ///Compute the state_age.
        ScalarType state_age = (ScalarType)i * dt;
        vec_contribution_to_foi_1[i] =
            parameters.get<TransitionDistributions>()[(int)InfectionTransition::ExposedToInfectedNoSymptoms].eval(
                state_age);
    }
    m_transitiondistributions[(int)InfectionState::Exposed] = vec_contribution_to_foi_1;

    // Other states:
    // Vector containing the transitions for the states with two outgoing flows.
    std::vector<std::vector<int>> vector_transitions = {
        {(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
         (int)InfectionTransition::InfectedNoSymptomsToRecovered},
        {(int)InfectionTransition::InfectedSymptomsToInfectedSevere,
         (int)InfectionTransition::InfectedSymptomsToRecovered},
        {(int)InfectionTransition::InfectedSevereToInfectedCritical,
         (int)InfectionTransition::InfectedSevereToRecovered},
        {(int)InfectionTransition::InfectedCriticalToDead, (int)InfectionTransition::InfectedCriticalToRecovered}};

    for (int state = 2; state < (int)InfectionState::Count - 2; state++) {
        Eigen::Index calc_time_index =
            (Eigen::Index)std::ceil(std::max(m_transitiondistributions_support_max[vector_transitions[state - 2][0]],
                                             m_transitiondistributions_support_max[vector_transitions[state - 2][1]]) /
                                    dt);

        // Create vec_derivative that stores the value of the approximated derivative for all necessary time points.
        std::vector<ScalarType> vec_contribution_to_foi_2(calc_time_index + 1, 0.);
        for (Eigen::Index i = 0; i <= calc_time_index; i++) {
            ///Compute the state_age.
            ScalarType state_age = (ScalarType)i * dt;
            vec_contribution_to_foi_2[i] =
                parameters.get<TransitionProbabilities>()[vector_transitions[state - 2][0]] *
                    parameters.get<TransitionDistributions>()[vector_transitions[state - 2][0]].eval(state_age) +
                parameters.get<TransitionProbabilities>()[vector_transitions[state - 2][1]] *
                    parameters.get<TransitionDistributions>()[vector_transitions[state - 2][1]].eval(state_age);
        }
        m_transitiondistributions[state] = vec_contribution_to_foi_2;
    }
}

void Model::set_transitiondistributions_derivative(ScalarType dt)
{
    // The transition SusceptibleToExposed is not needed in the computations.
    for (int transition = 1; transition < (int)InfectionTransition::Count; transition++) {
        Eigen::Index support_max_index =
            (Eigen::Index)std::ceil(m_transitiondistributions_support_max[transition] / dt);

        // Create vec_derivative that stores the value of the approximated derivative for all necessary time points.
        std::vector<ScalarType> vec_derivative(support_max_index + 1, 0.);

        for (Eigen::Index i = 0; i <= support_max_index; i++) {
            // Compute state_age.
            ScalarType state_age = (ScalarType)i * dt;
            //Compute the apprximate derivative by a backwards difference scheme.
            vec_derivative[i] = (parameters.get<TransitionDistributions>()[transition].eval(state_age) -
                                 parameters.get<TransitionDistributions>()[transition].eval(state_age - dt)) /
                                dt;
        }
        m_transitiondistributions_derivative[transition] = vec_derivative;
    }
}

void Model::set_meaninfectivity(ScalarType dt)
{

    // Determine the corresponding time index.
    Eigen::Index calc_time = std::min(
        m_calctime, m_transitiondistributions_support_max[(int)InfectionTransition::ExposedToInfectedNoSymptoms]);
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / dt) - 1;

    m_meaninfectivity = std::vector<ScalarType>(calc_time_index + 1, 0.);

    // The mean infectivity is computed as exp(-NaturalDeathRate t) * (a_2(t)-a_1(t)).
    for (int time_point_index = 1; time_point_index <= calc_time_index; time_point_index++) {
        ScalarType a_1 = 0;
        ScalarType a_2 = 0;

        ScalarType time_point_age = (ScalarType)time_point_index * dt;
        //We compute a_1 and a_2 at time_point.
        for (int i = 0; i <= time_point_index - 1; i++) {
            ScalarType state_age_i = static_cast<ScalarType>(i) * dt;
            //Computation of a_1.
            a_1 += parameters.get<RelativeTransmissionNoSymptoms>().eval(state_age_i) *
                   m_transitiondistributions[static_cast<int>(InfectionState::InfectedNoSymptoms)][i] *
                   m_transitiondistributions_derivative[(int)InfectionTransition::ExposedToInfectedNoSymptoms]
                                                       [time_point_index - i];
            //Computation of a_2.
            ScalarType sum = 0;
            for (int j = 0; j <= i - 1; j++) {
                ScalarType state_age_j = static_cast<ScalarType>(j) * dt;
                sum +=
                    parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age_j) *
                    m_transitiondistributions[static_cast<int>(InfectionState::InfectedSymptoms)][j] *
                    parameters.get<TransitionProbabilities>()[(
                        int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
                    m_transitiondistributions_derivative[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
                                                        [time_point_index - j];
            }
            a_2 += m_transitiondistributions_derivative[(int)InfectionTransition::ExposedToInfectedNoSymptoms]
                                                       [time_point_index - i] *
                   sum;
        }
        m_meaninfectivity[time_point_index] =
            std::exp(-parameters.get<NaturalBirthRate>() * time_point_age) * (dt * dt * a_2 - dt * a_1);
    }
}
// void Model::set_initalvaluesforceofinfection(ScalarType dt)
// {
//     // Determine the relevant calculation area.
//     ScalarType calc_time =
//         std::max({m_transitiondistributions_support_max[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms],
//                   m_transitiondistributions_support_max[(int)InfectionTransition::InfectedNoSymptomsToRecovered],
//                   m_transitiondistributions_support_max[(int)InfectionTransition::InfectedSymptomsToInfectedSevere],
//                   m_transitiondistributions_support_max[(int)InfectionTransition::InfectedSymptomsToRecovered]});

//     Eigen::Index calc_time_index   = (Eigen::Index)std::ceil(calc_time / dt) - 1;
//     m_initalvaluesforceofinfection = std::vector<ScalarType>(calc_time_index, 0.);

//     for (int time_point_index = 1; time_point_index < calc_time_index; time_point_index++) {
//         ScalarType time_point_age = time_point_index / dt;
//         // Evtl wieder hier einfügen wenn ich weiß wie !
//         //We start by adding the phi_0 term.
//         // m_initalvaluesforceofinfection[time_point_index] =
//         //     m_meaninfectivity[time_point_index] *
//         //     (populations[static_cast<int>(InfectionState::InfectedNoSymptoms)][0] *
//         //          m_transitiondistributions_in_forceofinfection[0][time_point_index] +
//         //      populations[static_cast<int>(InfectionState::InfectedSymptoms)][0] *
//         //          m_transitiondistributions_in_forceofinfection[1][time_point_index]);

//         // Now we compute the f term, that is f=f_3-f_2-f_1.
//         ScalarType f_1 = 0.;
//         ScalarType f_2 = 0.;
//         ScalarType f_3 = 0.;
//         for (int i = 1; i <= time_point_index; i++) {
//             ScalarType state_age_i = static_cast<ScalarType>(i - 1) * dt;
//             //Compute sum for f_1:
//             f_1 +=
//                 m_transitiondistributions_derivative[static_cast<int>(InfectionTransition::ExposedToInfectedNoSymptoms)]
//                                                     [time_point_index + 1 - i] *
//                 parameters.get<RelativeTransmissionNoSymptoms>().eval(state_age_i) *
//                 m_transitiondistributions_in_forceofinfection[0][i - 1];

//             //Compute sum for f_2:
//             f_2 += m_transitiondistributions_derivative[static_cast<int>(
//                        InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)][time_point_index + 1 - i] *
//                    parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age_i) *
//                    m_transitiondistributions_in_forceofinfection[1][i - 1];

//             //Compute sum for f_3:
//             ScalarType sum = 0.;
//             for (int j = 1; j <= i; j++) {
//                 ScalarType state_age_j = static_cast<ScalarType>(j - 1) * dt;
//                 sum += parameters.get<TransitionProbabilities>()[static_cast<int>(
//                            InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] *
//                        m_transitiondistributions_derivative[static_cast<int>(
//                            InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)][time_point_index + 1 - j] *
//                        parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age_j) *
//                        m_transitiondistributions_in_forceofinfection[1][j - 1];
//             }
//             f_3 +=
//                 m_transitiondistributions_derivative[static_cast<int>(InfectionTransition::ExposedToInfectedNoSymptoms)]
//                                                     [time_point_index + 1 - i] *
//                 sum;
//         }
//         m_initalvaluesforceofinfection[time_point_index] +=
//             std::pow(dt, 2) * populations[static_cast<int>(InfectionState::Exposed)][0] *
//                 std::exp(-parameters.get<NaturalBirthRate>() * time_point_age) * f_3 -
//             dt * populations[static_cast<int>(InfectionState::InfectedNoSymptoms)][0] *
//                 std::exp(-parameters.get<NaturalBirthRate>() * time_point_age) * f_2 -
//             dt * populations[static_cast<int>(InfectionState::Exposed)][0] *
//                 std::exp(-parameters.get<NaturalBirthRate>() * time_point_age) * f_1;
//     }
// }

// ---- Functionality for the model analysis. ----

void Model::set_reproductionnumber_c(ScalarType dt)
{
    // Determine the corresponding time index.
    Eigen::Index calc_time = std::min(
        m_calctime, m_transitiondistributions_support_max[(int)InfectionTransition::ExposedToInfectedNoSymptoms]);
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / dt) - 1;
    ScalarType sum               = 0;
    for (int i = 0; i <= calc_time_index; i++) {
        sum += m_meaninfectivity[i];
    }
    m_reproductionnumber_c = parameters.get<TransmissionProbabilityOnContact>().eval(0) *
                             parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(0)(0, 0) * dt * sum;
}

void Model::set_probability_of_transition(ScalarType dt)
{

    for (int transition = 1; transition < (int)InfectionTransition::Count; transition++) {
        Eigen::Index support_max_index =
            (Eigen::Index)std::ceil(m_transitiondistributions_support_max[transition] / dt);
        ScalarType sum = 0.;
        for (int i = 0; i <= support_max_index; i++) {
            ScalarType state_age_i = static_cast<ScalarType>(i) * dt;
            sum -= std::exp(-parameters.get<NaturalDeathRate>() * state_age_i) *
                   parameters.get<TransitionProbabilities>()[transition] *
                   m_transitiondistributions_derivative[transition][i];
        }
        m_probabilityoftransition[transition] = dt * sum;
    }
}

//TODO schönere implementierung überlegen, die evtl m_transitiondistributions_in_forceofinfection verwendet!!!!
void Model::set_meansojourntime(ScalarType dt)
{
    std::vector<int> relevant_states = {(int)InfectionState::Exposed, (int)InfectionState::InfectedNoSymptoms,
                                        (int)InfectionState::InfectedSymptoms, (int)InfectionState::InfectedSevere,
                                        (int)InfectionState::InfectedCritical};

    ScalarType sum = 0;
    for (int contribution = 0; contribution < 5; contribution++) {
        int calc_time_index = m_transitiondistributions[relevant_states[contribution]].size();
        for (int i = 0; i < calc_time_index; i++) {
            ScalarType state_age_i = (ScalarType)i * dt;
            sum += m_transitiondistributions[relevant_states[contribution]][i] *
                   std::exp(-parameters.get<NaturalDeathRate>() * state_age_i);
        }

        m_meansojourntime[relevant_states[contribution]] = dt * sum;
    }
}

void Model::set_equilibrium()
{
    // Case of the disease free equilibrium.
    if (m_reproductionnumber_c < 1) {
        m_equilibriumforceofinfection                                                = 0;
        m_equilibriumnormalizedcompartments[(int)InfectionState::Susceptible]        = 1;
        m_equilibriumnormalizedcompartments[(int)InfectionState::Exposed]            = 0;
        m_equilibriumnormalizedcompartments[(int)InfectionState::InfectedNoSymptoms] = 0;
        m_equilibriumnormalizedcompartments[(int)InfectionState::InfectedSymptoms]   = 0;
        m_equilibriumnormalizedcompartments[(int)InfectionState::InfectedSevere]     = 0;
        m_equilibriumnormalizedcompartments[(int)InfectionState::InfectedCritical]   = 0;
        m_equilibriumnormalizedcompartments[(int)InfectionState::Recovered]          = 0;
    }
    // Case of the endemic equilibrium.

    if (m_reproductionnumber_c > 1) {

        //Set probabilities that we transition to a specific compartment during the disease.
        ScalarType transition_to_INS = m_probabilityoftransition[(int)InfectionTransition::ExposedToInfectedNoSymptoms];
        ScalarType transition_to_ISy =
            transition_to_INS *
            m_probabilityoftransition[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms];
        ScalarType transition_to_ISe =
            transition_to_ISy * m_probabilityoftransition[(int)InfectionTransition::InfectedSymptomsToInfectedSevere];
        ScalarType transition_to_IC =
            transition_to_ISe * m_probabilityoftransition[(int)InfectionTransition::InfectedSevereToInfectedCritical];
        ScalarType transition_to_D =
            transition_to_IC * m_probabilityoftransition[(int)InfectionTransition::InfectedCriticalToDead];
        ScalarType transition_to_R =
            m_probabilityoftransition[(int)InfectionTransition::InfectedNoSymptomsToRecovered] * transition_to_INS +
            m_probabilityoftransition[(int)InfectionTransition::InfectedSymptomsToRecovered] * transition_to_ISy +
            m_probabilityoftransition[(int)InfectionTransition::InfectedSevereToRecovered] * transition_to_ISe +
            m_probabilityoftransition[(int)InfectionTransition::InfectedCriticalToRecovered] * transition_to_IC;

        m_equilibriumforceofinfection =
            parameters.get<NaturalDeathRate>() * (m_reproductionnumber_c - 1) +
            (parameters.get<NaturalDeathRate>() * transition_to_D * (m_reproductionnumber_c - 1)) /
                (1 - transition_to_D);
        m_equilibriumnormalizedcompartments[(int)InfectionState::Susceptible] = 1 / m_reproductionnumber_c;
        m_equilibriumnormalizedcompartments[(int)InfectionState::Exposed] =
            m_equilibriumforceofinfection * m_equilibriumnormalizedcompartments[(int)InfectionState::Susceptible] *
            m_meansojourntime[(int)InfectionState::Exposed];
        m_equilibriumnormalizedcompartments[(int)InfectionState::InfectedNoSymptoms] =
            m_equilibriumforceofinfection * m_equilibriumnormalizedcompartments[(int)InfectionState::Susceptible] *
            transition_to_INS * m_meansojourntime[(int)InfectionState::InfectedNoSymptoms];
        m_equilibriumnormalizedcompartments[(int)InfectionState::InfectedSymptoms] =
            m_equilibriumforceofinfection * m_equilibriumnormalizedcompartments[(int)InfectionState::Susceptible] *
            transition_to_ISy * m_meansojourntime[(int)InfectionState::InfectedSymptoms];
        m_equilibriumnormalizedcompartments[(int)InfectionState::InfectedSevere] =
            m_equilibriumforceofinfection * m_equilibriumnormalizedcompartments[(int)InfectionState::Susceptible] *
            transition_to_ISe * m_meansojourntime[(int)InfectionState::InfectedSevere];
        m_equilibriumnormalizedcompartments[(int)InfectionState::InfectedCritical] =
            m_equilibriumforceofinfection * m_equilibriumnormalizedcompartments[(int)InfectionState::Susceptible] *
            transition_to_IC * m_meansojourntime[(int)InfectionState::InfectedCritical];
        m_equilibriumnormalizedcompartments[(int)InfectionState::Recovered] =
            m_equilibriumforceofinfection * m_equilibriumnormalizedcompartments[(int)InfectionState::Susceptible] *
            transition_to_R / parameters.get<NaturalDeathRate>();
    }
}

}; // namespace endisecir

} // namespace mio
