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
#include <cstddef>
#include <vector>

namespace mio
{
namespace endisecir
{
Model::Model(TimeSeries<ScalarType>&& states_init)
    : parameters{Parameters()}
    , transitions{TimeSeries<ScalarType>(Eigen::Index(InfectionTransition::Count))}
    , populations{std::move(states_init)}
{
    // Set flows at start time t0.
    // As we assume that all individuals have infectio age 0 at time t0, the flows at t0 are set to 0.
    transitions.add_time_point(
        0, TimeSeries<ScalarType>::Vector::Constant(static_cast<size_t>(InfectionTransition::Count), 0.));

    // Set population size at start timt t0.
    ScalarType init_populationsize = std::accumulate(populations[0].begin(), populations[0].end(), 0);
    m_totalpopulation.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, init_populationsize));
    // Set the force of infection at start time t0.
    m_forceofinfection.add_time_point(0, TimeSeries<ScalarType>::Vector::Constant(1, 0));
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
             (1 - parameters.get<NaturalDeathRate>()) +
         m_totalpopulation[num_time_points - 2][0] * parameters.get<NaturalBirthRate>()) /
        (1 + dt * m_forceofinfection[num_time_points - 2][0]);
}

void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow,
                         Eigen::Index idx_CurrentCompartment, Eigen::Index current_time_index, ScalarType dt)
{
    ScalarType sum = 0;

    for (Eigen::Index i = 0; i < current_time_index; i++) {
        sum += transitions[i + 1][idx_IncomingFlow] *
               std::exp((-parameters.get<NaturalDeathRate>()) * (current_time_index - i)) *
               m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index - i];
    }
    transitions.get_value(current_time_index)[idx_InfectionTransitions] =
        (-dt) * parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum -
        std::exp((-parameters.get<NaturalDeathRate>()) * (current_time_index)) *
            populations[0][idx_CurrentCompartment] *
            parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] *
            m_transitiondistributions_derivative[idx_InfectionTransitions][current_time_index];
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
        dt * m_forceofinfection[current_time_index - 1][0] *
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

void Model::update_compartment_from_flow(InfectionState infectionState,
                                         std::vector<InfectionTransition> const& IncomingFlows,
                                         std::vector<InfectionTransition> const& OutgoingFlows,
                                         bool NaturalDeathispossible, ScalarType dt)
{
    Eigen::Index num_time_points = populations.get_num_time_points();

    ScalarType updated_compartment = 0;
    if (NaturalDeathispossible) {
        updated_compartment = populations[num_time_points - 2][static_cast<int>(infectionState)] *
                              (1 - parameters.get<NaturalDeathRate>());
    }
    else {
        updated_compartment = populations[num_time_points - 2][static_cast<int>(infectionState)];
    }
    for (const InfectionTransition& inflow : IncomingFlows) {
        updated_compartment += dt * transitions.get_last_value()[(int)inflow];
    }
    for (const InfectionTransition& outflow : OutgoingFlows) {
        updated_compartment -= dt * transitions.get_last_value()[(int)outflow];
    }
    populations.get_last_value()[(int)infectionState] = updated_compartment;
}

void Model::update_compartments(ScalarType dt)
{
    // Exposed
    update_compartment_from_flow(InfectionState::Exposed, {InfectionTransition::SusceptibleToExposed},
                                 {InfectionTransition::ExposedToInfectedNoSymptoms}, true, dt);

    // InfectedNoSymptoms
    update_compartment_from_flow(
        InfectionState::InfectedNoSymptoms, {InfectionTransition::ExposedToInfectedNoSymptoms},
        {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms, InfectionTransition::InfectedNoSymptomsToRecovered},
        true, dt);

    // InfectedSymptoms
    update_compartment_from_flow(
        InfectionState::InfectedSymptoms, {InfectionTransition::InfectedNoSymptomsToInfectedSymptoms},
        {InfectionTransition::InfectedSymptomsToInfectedSevere, InfectionTransition::InfectedSymptomsToRecovered}, true,
        dt);

    // InfectedSevere
    update_compartment_from_flow(
        InfectionState::InfectedSevere, {InfectionTransition::InfectedSymptomsToInfectedSevere},
        {InfectionTransition::InfectedSevereToInfectedCritical, InfectionTransition::InfectedSevereToRecovered}, true,
        dt);

    // InfectedCritical
    update_compartment_from_flow(
        InfectionState::InfectedCritical, {InfectionTransition::InfectedSevereToInfectedCritical},
        {InfectionTransition::InfectedCriticalToDead, InfectionTransition::InfectedCriticalToRecovered}, true, dt);

    // Recovered
    update_compartment_from_flow(InfectionState::Recovered,
                                 {
                                     InfectionTransition::InfectedNoSymptomsToRecovered,
                                     InfectionTransition::InfectedSymptomsToRecovered,
                                     InfectionTransition::InfectedSevereToRecovered,
                                     InfectionTransition::InfectedCriticalToRecovered,
                                 },
                                 std::vector<InfectionTransition>(), true, dt);

    // Dead
    update_compartment_from_flow(InfectionState::Dead, {InfectionTransition::InfectedCriticalToDead},
                                 std::vector<InfectionTransition>(), false, dt);
}

void Model::compute_populationsize()
{
    m_totalpopulation.get_last_value()[0] =
        std::accumulate(populations.get_last_value().begin(), populations.get_last_value().end(), 0);
}

void Model::compute_forceofinfection(ScalarType dt)
{

    Eigen::Index num_time_points = populations.get_num_time_points();
    ScalarType current_time      = populations.get_last_time();

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(m_calctime / dt) - 1;
    ScalarType sum               = 0;

    //Computation of phi_0
    ScalarType phi_0 = 0;
    if (num_time_points <= calc_time_index) {
        phi_0 = m_meaninfectivity[num_time_points] *
                (populations[0][static_cast<int>(InfectionState::InfectedNoSymptoms)] *
                     m_transitiondistributions_in_forceofinfection[0][num_time_points] +
                 populations[0][static_cast<int>(InfectionState::InfectedSymptoms)] *
                     m_transitiondistributions_in_forceofinfection[1][num_time_points]) /
                m_totalpopulation[num_time_points - 1][0];
    }
    //Computation of f
    ScalarType f   = 0.;
    ScalarType f_1 = 0.;
    ScalarType f_2 = 0.;
    ScalarType f_3 = 0.;
    for (int i = 1; i <= calc_time_index; i++) {
        ScalarType state_age_i = static_cast<ScalarType>(i - 1) * dt;
        //Compute sum for f_1:
        f_1 += m_transitiondistributions_derivative[static_cast<int>(InfectionTransition::ExposedToInfectedNoSymptoms)]
                                                   [current_time + 1 - i] *
               parameters.get<RelativeTransmissionNoSymptoms>().eval(state_age_i) *
               m_transitiondistributions_in_forceofinfection[0][i - 1];

        //Compute sum for f_2:
        f_2 += m_transitiondistributions_derivative[static_cast<int>(
                   InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)][num_time_points + 1 - i] *
               parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age_i) *
               m_transitiondistributions_in_forceofinfection[1][i - 1];

        //Compute sum for f_3:
        sum = 0.;
        for (int j = 1; j <= i; j++) {
            ScalarType state_age_j = static_cast<ScalarType>(j - 1) * dt;
            sum += parameters.get<TransitionProbabilities>()[static_cast<int>(
                       InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] *
                   m_transitiondistributions_derivative[static_cast<int>(
                       InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)][num_time_points + 1 - j] *
                   parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age_j) *
                   m_transitiondistributions_in_forceofinfection[1][j - 1];
        }
        f_3 += m_transitiondistributions_derivative[static_cast<int>(InfectionTransition::ExposedToInfectedNoSymptoms)]
                                                   [num_time_points + 1 - i] *
               sum;
    }
    f = std::pow(dt, 2) * populations[0][static_cast<int>(InfectionState::Exposed)] *
            std::exp(-parameters.get<NaturalBirthRate>() * current_time) * f_3 -
        dt * populations[0][static_cast<int>(InfectionState::InfectedNoSymptoms)] *
            std::exp(-parameters.get<NaturalBirthRate>() * current_time) * f_2 -
        dt * populations[0][static_cast<int>(InfectionState::Exposed)] *
            std::exp(-parameters.get<NaturalBirthRate>() * current_time) * f_1;

    // Computation of the rest

    //ERKLÄRUNG!!!
    Eigen::Index starting_point = std::max(0, (int)num_time_points - (int)calc_time_index - 1);

    sum = 0.;
    for (Eigen::Index i = starting_point; i < num_time_points - 1; i++) {
        std::cout << " mean infectivitity at time " << num_time_points - 1 - i << " is "
                  << m_meaninfectivity[num_time_points - 1 - i] << "\n";
        sum += m_meaninfectivity[num_time_points - 1 - i] * populations[i + 1][(int)InfectionState::Susceptible] *
               m_forceofinfection[i][0];
    }

    m_forceofinfection.get_last_value()[0] =
        phi_0 + (dt * sum + f) *
                    (parameters.get<TransmissionProbabilityOnContact>().eval(current_time) *
                     parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(num_time_points)(0, 0)) /
                    m_totalpopulation[num_time_points - 1][0];
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

    // Transitions needed in the force of infection term.
    std::vector<std::vector<int>> relevant_transitions = {
        {(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
         (int)InfectionTransition::InfectedNoSymptomsToRecovered},
        {(int)InfectionTransition::InfectedSymptomsToInfectedSevere,
         (int)InfectionTransition::InfectedSymptomsToRecovered}};

    // Determine the corresponding time index to m_calctime.
    // Subtract 1 because in the last summand all TransitionDistributions evaluate to 0 (by definition of support_max).
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(m_calctime / dt) - 1;

    // Compute the distributions from survival functions and transition probabilities starting in InfectedNoSymptoms and
    // InfectedSymptims ont the force of infection term.
    for (int contribution = 0; contribution < 2; contribution++) {
        std::vector<ScalarType> vec_contribution_to_foi(calc_time_index + 1, 0.);
        for (Eigen::Index i = 0; i <= calc_time_index; i++) {
            //Compute the state_age.
            ScalarType state_age = i * dt;
            vec_contribution_to_foi[i] =
                parameters.get<TransitionProbabilities>()[relevant_transitions[contribution][0]] *
                    parameters.get<TransitionDistributions>()[relevant_transitions[contribution][0]].eval(state_age) +
                parameters.get<TransitionProbabilities>()[relevant_transitions[contribution][1]] *
                    parameters.get<TransitionDistributions>()[relevant_transitions[contribution][1]].eval(state_age);
        }
        m_transitiondistributions_in_forceofinfection[contribution] = vec_contribution_to_foi;
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
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(m_calctime / dt) - 1;

    m_meaninfectivity = std::vector<ScalarType>(calc_time_index, 0.);

    // The mean infectivity is computed as exp(-NaturalDeathRate t) * (a_2(t)-a_1(t)).
    for (int time_point_index = 1; time_point_index <= calc_time_index; time_point_index++) {
        ScalarType a_1 = 0;
        ScalarType a_2 = 0;

        ScalarType time_point_age = (ScalarType)time_point_index * dt;
        //We compute a_1 and a_2 at time_point.
        for (int i = 1; i <= time_point_index; i++) {
            ScalarType state_age_i = static_cast<ScalarType>(i) * dt;
            //Computation of a_1.
            a_1 += parameters.get<RelativeTransmissionNoSymptoms>().eval(state_age_i) *
                   m_transitiondistributions_in_forceofinfection[0][i] *
                   m_transitiondistributions_derivative[(int)InfectionTransition::ExposedToInfectedNoSymptoms]
                                                       [time_point_index - i];
            //Computation of a_2.
            ScalarType sum = 0;
            for (int j = 1; j <= i; j++) {
                ScalarType state_age_j = static_cast<ScalarType>(j) * dt;
                sum +=
                    parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age_j) *
                    m_transitiondistributions_in_forceofinfection[1][j] *
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
            std::exp(-parameters.get<NaturalBirthRate>() * time_point_age) * (std::pow(dt, 2) * a_2 - dt * a_1);
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

// ---- Computation of the Reproduction number. ----

void Model::set_reproductionnumber_c(ScalarType dt)
{
    ScalarType sum = 0;
    for (ScalarType mean_infectivity_i : m_meaninfectivity) {
        sum += mean_infectivity_i;
    }
    m_reproductionnumber_c = parameters.get<TransmissionProbabilityOnContact>().eval(0) *
                             parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(0)(0, 0) * dt * sum;
}

}; // namespace endisecir

} // namespace mio
