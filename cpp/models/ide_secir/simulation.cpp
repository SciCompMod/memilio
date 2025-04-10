/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Martin J Kuehn, Anna Wendler, Lena Ploetzke
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
#include "ide_secir/simulation.h"
#include "ide_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include <memory>

#include <iostream>
#include <cmath>
#include <chrono>
#include <openacc.h>


namespace mio
{
namespace isecir
{

void Simulation::advance(ScalarType tmax)
{
    mio::log_info("Simulating IDE-SECIR from t0 = {} until tmax = {} with dt = {}.",
                  m_model->transitions.get_last_time(), tmax, m_dt);
    m_model->set_transitiondistributions_support_max(m_dt);
    m_model->set_transitiondistributions_derivative(m_dt);
    m_model->set_transitiondistributions_in_forceofinfection(m_dt);
    m_model->initial_compute_compartments(m_dt);


    // For every time step:
    while (m_model->transitions.get_last_time() < tmax - m_dt / 2) {

        m_model->transitions.add_time_point(m_model->transitions.get_last_time() + m_dt);
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt);

        // compute Susceptibles:
        m_model->compute_susceptibles(m_dt);

        // Compute flows:
        m_model->flows_current_timestep(m_dt);

        // Update remaining compartments:
        m_model->update_compartments();

        // // Compute m_forceofinfection (only used for calculation of Susceptibles and flow SusceptibleToExposed in the next timestep!):
        // m_model->compute_forceofinfection(m_dt);

        // Measure ONLY the time of compute_forceofinfection()
        auto start = std::chrono::high_resolution_clock::now(); // Start timer
        m_model->compute_forceofinfection(m_dt);                // Function to measure
        auto end = std::chrono::high_resolution_clock::now();   // Stop timer

        // Calculate duration in microseconds (or milliseconds if preferred)
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        // Print ONLY the time taken by compute_forceofinfection()
        std::cout << "compute_forceofinfection() took: " << duration.count() << " μs" << std::endl;
    }
}











































// namespace mio
// {
// namespace isecir
// {

// void Simulation::advance(ScalarType tmax)
// {
//     mio::log_info("Simulating IDE-SECIR from t0 = {} until tmax = {} with dt = {}.",
//                   m_model->transitions.get_last_time(), tmax, m_dt);
//     m_model->set_transitiondistributions_support_max(m_dt);
//     m_model->set_transitiondistributions_derivative(m_dt);
//     m_model->set_transitiondistributions_in_forceofinfection(m_dt);
//     m_model->initial_compute_compartments(m_dt);


//     TimeSeries<ScalarType>& trans = m_model->transitions;
//     TimeSeries<ScalarType>& pop = m_model->populations;
//     Parameters& par = m_model->parameters;

//     int number_transitions = (int)mio::isecir::InfectionTransition::Count;
//     int number_compartments = (int)mio::isecir::InfectionState::Count;

//     // CustomIndexArray<ScalarType, AgeGroup>& forceofinfection = m_model->m_forceofinfection;
//     int number_age_groups = m_model->m_num_agegroups;
//     std::vector<ScalarType> v_forceofinfection(number_age_groups, 0.0);
//     for (int i = 0; i < number_age_groups; i++){
//         v_forceofinfection[i] = m_model->m_forceofinfection.begin()[i];
//     }
    
//     std::vector<std::vector<ScalarType>> v_trans(
//         trans.get_num_time_points(), 
//         std::vector<ScalarType>(trans.get_num_elements())
//     );
//     for (int i = 0; i < trans.get_num_time_points(); i++) {
//         for (int j = 0; j < trans.get_num_elements(); j++) {
//             v_trans[i][j] = trans[i][j];
//             // std::cout<< i << ", " << j << ": "<< v_trans[i][j]<<std::endl;
//         }
//     }


//     std::vector<std::vector<ScalarType>> v_pop(
//         pop.get_num_time_points(), 
//         std::vector<ScalarType>(pop.get_num_elements())
//     );
//     for (int i = 0; i < pop.get_num_time_points(); i++) {
//         for (int j = 0; j < pop.get_num_elements(); j++) {
//             v_pop[i][j] = pop[i][j];
//             std::cout<< i << ", " << j << ": "<< v_pop[i][j]<<std::endl;
//         }
//     }



 
//     // Begin While Loop!


//     // m_model->transitions.add_time_point(m_model->transitions.get_last_time() + m_dt);
//     std::vector<ScalarType> zero_trans(v_trans[0].size(), 0.0);
//     v_trans.push_back(zero_trans);
//     // m_model->transitions.add_time_point(m_model->transitions.get_last_time() + m_dt);
//     std::vector<ScalarType> zero_pop(v_pop[0].size(), 0.0); 
//     v_pop.push_back(zero_pop);



//     // Funktion: m_model->compute_susceptibles(m_dt):
//     // for (int i = 0; i < (m_num_agegroups); ++i) {
//     //     AgeGroup group(i);
//     //     Eigen::Index num_time_points = populations.get_num_time_points();
//     //     int Si                           = get_state_flat_index(Eigen::Index(InfectionState::Susceptible), group);
//     //     populations.get_last_value()[Si] = populations[num_time_points - 2][Si] / (1 + dt * m_forceofinfection[group]);
//     // }


//     // std::cout<<"Start: "<<std::endl;
//     size_t prev_time_idx = v_pop.size() - 2; 
//     int num_time_points = v_pop.size();
//     for (int i = 0; i < number_age_groups; ++i) {
//         int Si = i*number_compartments;
//         // Update susceptible compartment in latest time point
//         v_pop.back()[Si] = v_pop[prev_time_idx][Si] / (1 + m_dt * v_forceofinfection[i]);
//         // std::cout<<v_pop.back()[Si]<<std::endl;
//     }


//     // Funktion: m_model->flows_current_timestep(m_dt):
//     // for (AgeGroup group = AgeGroup(0); group < AgeGroup(m_num_agegroups); ++group) {
//     //     int StEi = get_transition_flat_index(Eigen::Index(InfectionTransition::SusceptibleToExposed), group);
//     //     int Si   = get_state_flat_index(Eigen::Index(InfectionState::Susceptible), group);

//     //     // Calculate flow SusceptibleToExposed with force of infection from previous time step and Susceptibles from current time step.
//     //     transitions.get_last_value()[StEi] = dt * m_forceofinfection[group] * populations.get_last_value()[Si];

//     //     // Calculate all other flows with compute_flow.
//     //     // Flow from Exposed to InfectedNoSymptoms
//     //     compute_flow(Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
//     //                  Eigen::Index(InfectionTransition::SusceptibleToExposed), dt, group);
//     //     // Flow from InfectedNoSymptoms to InfectedSymptoms
//     //     compute_flow(Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
//     //                  Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, group);
//     //     // Flow from InfectedNoSymptoms to Recovered
//     //     compute_flow(Eigen::Index(InfectionTransition::InfectedNoSymptomsToRecovered),
//     //                  Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt, group);
//     //     // Flow from InfectedSymptoms to InfectedSevere
//     //     compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
//     //                  Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt, group);
//     //     // Flow from InfectedSymptoms to Recovered
//     //     compute_flow(Eigen::Index(InfectionTransition::InfectedSymptomsToRecovered),
//     //                  Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt, group);
//     //     // Flow from InfectedSevere to InfectedCritical
//     //     compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
//     //                  Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, group);
//     //     // Flow from InfectedSevere to Recovered
//     //     compute_flow(Eigen::Index(InfectionTransition::InfectedSevereToRecovered),
//     //                  Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt, group);
//     //     // Flow from InfectedCritical to Dead
//     //     compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToDead),
//     //                  Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, group);
//     //     // Flow from InfectedCritical to Recovered
//     //     compute_flow(Eigen::Index(InfectionTransition::InfectedCriticalToRecovered),
//     //                  Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt, group);
//     // }

//     //     void Model::compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt,
//     //         Eigen::Index current_time_index, AgeGroup group)
//     // {
//     // ScalarType sum = 0;
//     // /* In order to satisfy TransitionDistribution(dt*i) = 0 for all i >= k, k is determined by the maximum support of 
//     // the distribution.
//     // Since we are using a backwards difference scheme to compute the derivative, we have that the
//     // derivative of TransitionDistribution(dt*i) = 0 for all i >= k+1.

//     // Hence calc_time_index goes until std::ceil(support_max/dt) since for std::ceil(support_max/dt)+1 all terms are 
//     // already zero. 
//     // This needs to be adjusted if we are changing the finite difference scheme */
//     // Eigen::Index calc_time_index =
//     // (Eigen::Index)std::ceil(m_transitiondistributions_support_max[group][idx_InfectionTransitions] / dt);

//     // int transition_idx = get_transition_flat_index(idx_InfectionTransitions, size_t(group));
//     // for (Eigen::Index i = current_time_index - calc_time_index; i < current_time_index; i++) {
//     // // (current_time_index - i)  is the index corresponding to time the individuals have already spent in this state.
//     // sum += m_transitiondistributions_derivative[group][idx_InfectionTransitions][current_time_index - i] *
//     // transitions[i + 1][idx_IncomingFlow];
//     // }

//     // transitions.get_value(current_time_index)[transition_idx] =
//     //     (-dt) * parameters.get<TransitionProbabilities>()[group][idx_InfectionTransitions] * sum;
//     // }









//     // Print v_pop back








//     // % Lese aus par die Werte in die vektoren ein

//     // TransitionDistributions
//     // TransitionProbabilities (vec<vec<vec<float>>>)
//     // ContactPatterns (vec<vec<float>>)
//     // TransmissionProbabilityOnContact (vec<float>)
//     // RelativeTransmissionNoSymptoms (vec<float>)
//     // RiskOfInfectionFromSymptomatic (vec<float>)
//     // StartDay (float)
//     // Seasonality (float)




















//     // // ScalarType* par_ptr = ...
//     // // std::size_t par_size = ...


//     // // ScalarType* trans_ptr = trans.data();
//     // // ScalarType* pop_ptr = pop.data();
//     // // ScalarType* trans_ptr = trans.data();

//     // // std::cout<<trans_ptr[0]<<std::endl;




//     // // For every time step:
//     // while (m_model->transitions.get_last_time() < tmax - m_dt / 2) {

//     //     m_model->transitions.add_time_point(m_model->transitions.get_last_time() + m_dt);
//     //     m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt);

//     //     // compute Susceptibles:
//     //     m_model->compute_susceptibles(m_dt);

//     //     // Compute flows:
//     //     m_model->flows_current_timestep(m_dt);

//     //     // Update remaining compartments:
//     //     m_model->update_compartments();

//     //     // Compute m_forceofinfection (only used for calculation of Susceptibles and flow SusceptibleToExposed in the next timestep!):
//     //     m_model->compute_forceofinfection(m_dt);
//     // }
// }

TimeSeries<ScalarType> simulate(ScalarType tmax, ScalarType dt, Model const& m_model)
{
    m_model.check_constraints(dt);
    Simulation sim(m_model, dt);
    sim.advance(tmax);
    return sim.get_result();
}

} // namespace isecir
} // namespace mio
