#include "ide_endemic_secir/normalized_model.h"
#include "ide_endemic_secir/computed_parameters.h"
#include "ide_endemic_secir/model.h"
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
#include <cstdlib>
#include <memory>
#include <vector>

namespace mio
{
namespace endisecir
{

NormModel::NormModel(CompParameters const& compparams)
    : compparameters{std::make_shared<CompParameters>(compparams)}
    , transitions{TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count))}
    , populations{TimeSeries<ScalarType>(Eigen::Index(InfectionState::Count) - 1)}
{
    //Set populations at start time t0.
    TimeSeries<ScalarType>::Vector vec_initnormalizedpopulations =
        TimeSeries<ScalarType>::Vector(Eigen::Index(InfectionState::Count) - 1);
    for (int infection_state = 0; infection_state < Eigen::Index(InfectionState::Count) - 1; infection_state++) {
        vec_initnormalizedpopulations[infection_state] =
            compparameters->m_statesinit[0][infection_state] / compparameters->m_totalpopulationinit;
    }
    populations.add_time_point(0, vec_initnormalizedpopulations);

    // Set flows at start time t0.
    // As we assume that all individuals have infectio age 0 at time t0, the flows at t0 are set to 0.
    transitions.add_time_point(
        0, TimeSeries<ScalarType>::Vector::Constant(static_cast<size_t>(InfectionTransition::Count), 0.));

    m_forceofinfection.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, compparameters->m_NormFoI_0[0]));
}

bool NormModel::check_constraints() const
{
    if (!(static_cast<size_t>(populations.get_num_elements()) == static_cast<size_t>(InfectionState::Count) - 1)) {
        log_error(" A variable given for model construction is not valid. Number of elements in vector of populations "
                  "does not match the required number.");
        return true;
    }

    for (int i = 0; i < static_cast<int>(InfectionState::Count) - 1; i++) {
        if (populations[0][i] < 0) {
            log_error("Initialization failed. Initial values for populations are less than zero.");
            return true;
        }
    }
    return compparameters->check_constraints();
}
// ----Functionality for Initialization. ----

// ----Functionality for the iterations of a simulation. ----
void NormModel::compute_susceptibles(ScalarType dt)
{
    Eigen::Index num_time_points = populations.get_num_time_points();
    populations.get_last_value()[static_cast<int>(InfectionState::Susceptible)] =
        (populations[num_time_points - 2][static_cast<int>(InfectionState::Susceptible)] +
         dt * compparameters->parameters.get<NaturalBirthRate>()) /
        (1 + dt * (m_forceofinfection[num_time_points - 2][0] + compparameters->parameters.get<NaturalBirthRate>()) -
         transitions[num_time_points - 2][(Eigen::Index)InfectionTransition::InfectedCriticalToDead]);
}

void NormModel::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow,
                             Eigen::Index idx_CurrentCompartment, Eigen::Index current_time_index, ScalarType dt)
{
    Eigen::Index calc_time_index =
        (Eigen::Index)std::ceil(compparameters->m_transitiondistributions_support_max[idx_InfectionTransitions] / dt);
    ScalarType current_time_age = static_cast<ScalarType>(current_time_index) * dt;
    ScalarType sum              = 0;
    //Determine the starting point of the for loop.
    Eigen::Index starting_point = std::max(0, (int)current_time_index - (int)calc_time_index);

    for (Eigen::Index i = starting_point; i < current_time_index; i++) {
        ScalarType state_age_i = static_cast<ScalarType>(i) * dt;

        sum += (transitions[i + 1][idx_IncomingFlow] +
                (compparameters->parameters.get<NaturalDeathRate>() +
                 transitions[i][(Eigen::Index)InfectionTransition::InfectedCriticalToDead] -
                 compparameters->parameters.get<NaturalBirthRate>()) *
                    populations[i][idx_CurrentCompartment]) *
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

void NormModel::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow,
                             Eigen::Index idx_CurrentCompartment, ScalarType dt)
{
    Eigen::Index current_time_index = transitions.get_num_time_points() - 1;
    compute_flow(idx_InfectionTransitions, idx_IncomingFlow, idx_CurrentCompartment, current_time_index, dt);
}

void NormModel::flows_currents_timestep(ScalarType dt)
{
    Eigen::Index current_time_index = populations.get_num_time_points() - 1;
    // Calculate the transition SusceptibleToExposed with force of infection from previous time step and Susceptibles from
    // current time step.
    transitions.get_last_value()[static_cast<int>(InfectionTransition::SusceptibleToExposed)] =
        m_forceofinfection[current_time_index - 1][0] *
        populations.get_last_value()[static_cast<int>(InfectionState::Susceptible)];

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

void NormModel::update_compartment_with_sum(InfectionState infectionState,
                                            std::vector<InfectionTransition> const& IncomingFlows,
                                            bool Transitionispossible, ScalarType dt)
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
        if (Transitionispossible) {
            sum += compparameters->m_transitiondistributions[(int)infectionState][current_time_index - i] *
                   (sum_inflows + (compparameters->parameters.get<NaturalDeathRate>() +
                                   transitions[i][(Eigen::Index)InfectionTransition::InfectedCriticalToDead] -
                                   compparameters->parameters.get<NaturalBirthRate>()) *
                                      populations[i][(int)infectionState]) *
                   std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i));
        }
        else {
            sum += (sum_inflows + (compparameters->parameters.get<NaturalDeathRate>() +
                                   transitions[i][(Eigen::Index)InfectionTransition::InfectedCriticalToDead] -
                                   compparameters->parameters.get<NaturalBirthRate>()) *
                                      populations[i][(int)infectionState]) *
                   std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age - state_age_i));
        }
    }
    if (Transitionispossible) {
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
    else {
        populations.get_last_value()[(int)infectionState] =
            dt * sum + populations[0][(int)infectionState] *
                           std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time_age));
    }
}

void NormModel::update_compartments(ScalarType dt)
{

    // Exposed
    update_compartment_with_sum(InfectionState::Exposed, {InfectionTransition::SusceptibleToExposed}, true, dt);
    // InfectedNoSymptoms
    update_compartment_with_sum(InfectionState::InfectedNoSymptoms, {InfectionTransition::ExposedToInfectedNoSymptoms},
                                true, dt);

    // InfectedSymptoms
    update_compartment_with_sum(InfectionState::InfectedSymptoms,
                                {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms}, true, dt);

    // InfectedSevere
    update_compartment_with_sum(InfectionState::InfectedSevere, {InfectionTransition::InfectedSymptomsToInfectedSevere},
                                true, dt);

    // InfectedCritical
    update_compartment_with_sum(InfectionState::InfectedCritical,
                                {InfectionTransition::InfectedSevereToInfectedCritical}, true, dt);
    // Recovered
    update_compartment_with_sum(InfectionState::Recovered,
                                {
                                    InfectionTransition::InfectedNoSymptomsToRecovered,
                                    InfectionTransition::InfectedSymptomsToRecovered,
                                    InfectionTransition::InfectedSevereToRecovered,
                                    InfectionTransition::InfectedCriticalToRecovered,
                                },
                                false, dt);
}

void NormModel::compute_forceofinfection(ScalarType dt)
{

    Eigen::Index num_time_points = populations.get_num_time_points();
    ScalarType current_time      = populations.get_last_time();

    // Determine the starting point of the for loop.
    Eigen::Index starting_point1 = std::max(0, (int)num_time_points - 1 - (int)compparameters->m_infectivity.size());

    ScalarType sum1 = 0.;
    // Compute the first sum in the force of infection term.
    for (Eigen::Index i = starting_point1; i < num_time_points - 1; i++) {
        sum1 += compparameters->m_infectivity[num_time_points - 1 - i] *
                (populations[i + 1][(int)InfectionState::Susceptible] * m_forceofinfection[i][0] +
                 (compparameters->parameters.get<NaturalDeathRate>() +
                  transitions[i + 1][(int)InfectionTransition::InfectedCriticalToDead] -
                  compparameters->parameters.get<NaturalBirthRate>()) *
                     populations[i + 1][(int)InfectionState::Exposed]);
    }

    // Determine the starting point of the for loop.
    Eigen::Index starting_point2 = std::max(0, (int)num_time_points - 1 - (int)compparameters->m_B.size());

    ScalarType sum2 = 0.;
    // Compute the second sum in the force of infection term.
    for (Eigen::Index i = starting_point2; i < num_time_points - 1; i++) {
        ScalarType state_age = static_cast<ScalarType>(i) * dt;
        sum2 -= (compparameters->parameters.get<NaturalDeathRate>() +
                 transitions[i + 1][(int)InfectionTransition::InfectedCriticalToDead] -
                 compparameters->parameters.get<NaturalBirthRate>()) *
                std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time - state_age)) *
                populations[i + 1][(int)InfectionState::InfectedNoSymptoms] *
                compparameters->m_B[num_time_points - 1 - i];
    }

    // Determine the starting point of the for loop.
    Eigen::Index starting_point3 =
        std::max(0, (int)num_time_points - 1 -
                        (int)compparameters->m_transitiondistributions[(int)InfectionState::InfectedNoSymptoms].size());

    ScalarType sum3 = 0.;
    // Compute the third sum in the force of infection term.
    for (Eigen::Index i = starting_point3; i < num_time_points - 1; i++) {
        ScalarType state_age = static_cast<ScalarType>(i) * dt;
        sum3 +=
            (compparameters->parameters.get<NaturalDeathRate>() +
             transitions[i + 1][(int)InfectionTransition::InfectedCriticalToDead] -
             compparameters->parameters.get<NaturalBirthRate>()) *
            std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time - state_age)) *
            populations[i + 1][(int)InfectionState::InfectedNoSymptoms] *
            compparameters->parameters.get<RelativeTransmissionNoSymptoms>().eval(current_time - state_age) *
            compparameters->m_transitiondistributions[(int)InfectionState::InfectedNoSymptoms][num_time_points - 1 - i];
    }

    // Determine the starting point of the for loop.
    Eigen::Index starting_point4 =
        std::max(0, (int)num_time_points - 1 -
                        (int)compparameters->m_transitiondistributions[(int)InfectionState::InfectedSymptoms].size());

    ScalarType sum4 = 0.;
    // Compute the third sum in the force of infection term.
    for (Eigen::Index i = starting_point4; i < num_time_points - 1; i++) {
        ScalarType state_age = static_cast<ScalarType>(i) * dt;
        sum4 +=
            (compparameters->parameters.get<NaturalDeathRate>() +
             transitions[i + 1][(int)InfectionTransition::InfectedCriticalToDead] -
             compparameters->parameters.get<NaturalBirthRate>()) *
            std::exp(-compparameters->parameters.get<NaturalDeathRate>() * (current_time - state_age)) *
            populations[i + 1][(int)InfectionState::InfectedSymptoms] *
            compparameters->parameters.get<RiskOfInfectionFromSymptomatic>().eval(current_time - state_age) *
            compparameters->m_transitiondistributions[(int)InfectionState::InfectedSymptoms][num_time_points - 1 - i];
    }
    ScalarType sum = sum1 + sum2 + sum3 + sum4;

    m_forceofinfection.get_last_value()[0] =
        dt * sum *
        (compparameters->parameters.get<TransmissionProbabilityOnContact>().eval(current_time) *
         compparameters->parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(0, 0));

    // Add inital functions for the force of infection in case they still have an influence.
    if (num_time_points <= (int)compparameters->m_NormFoI_0.size()) {
        m_forceofinfection.get_last_value()[0] += compparameters->m_NormFoI_0[num_time_points - 1];
    }
    if (num_time_points <= (int)compparameters->m_NormInitFoI.size()) {
        m_forceofinfection.get_last_value()[0] += compparameters->m_NormInitFoI[num_time_points - 1];
    }
}

} // namespace endisecir

} // namespace mio