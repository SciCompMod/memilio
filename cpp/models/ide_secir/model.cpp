/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler, Lena Ploetzke
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
#include "infection_state.h"
#include <cstdio>
#include <iostream>

namespace mio
{
namespace isecir
{

Model::Model(TimeSeries<ScalarType>&& init, ScalarType dt_init, size_t N_init, size_t Dead0, Pa Parameterset_init)
    : parameters{Parameterset_init}
    , m_transitions{std::move(init)}
    , m_SECIR{TimeSeries<ScalarType>((double)InfectionState::Count)}
    , m_dt{dt_init}
    , m_N{N_init}
{
    m_SECIR.add_time_point(0);
    m_SECIR[Eigen::Index(0)][Eigen::Index(InfectionState::Dead)] = Dead0;
}

// define helper functions later needed for simulation
void Model::update_susceptibles()
{   
    m_SECIR.add_time_point(m_SECIR.get_last_time() + m_dt);
    Eigen::Index num_time_points = m_SECIR.get_num_time_points();
    // using number of susceptibles from previous time step and force of infection from previous time step: compute current number of susceptibles and store susceptibles in m_SECIR
    m_SECIR[num_time_points - 1][Eigen::Index(InfectionState::Susceptible)]=
            m_SECIR[num_time_points - 2][Eigen::Index(InfectionState::Susceptible)] / (1 + m_dt * m_forceofinfection);
}

void Model::update_forceofinfection()
{
    // determine force of infection for the current last point of time in m_transitions, should also be the last
    m_forceofinfection = 0;
    ScalarType calc_time      = 0;
    // determine the relevant calculation area; just in the support of one of the transition distributions
    // kann man hier nicht auch eine Funktion nehmen, die das Maximum berechnet für eine Menge von Indizes? Gibt es sowas schon?
    for (int i = 1; i < 5; i++) {
        if (parameters.get<TransitionDistributions>()[i].get_xright() > calc_time) {
            calc_time = parameters.get<TransitionDistributions>()[i].get_xright();
        }
    }
    
    //determine corresponding indices
    // TODO: Parameter müssten zum teil noch abhängig von tau und t sein, das ist in parameters aber noch nicht
    Eigen::Index calc_time_index                 = (Eigen::Index)std::ceil(calc_time / m_dt); //need calc_time timesteps in sum

    Eigen::Index num_time_points = m_transitions.get_num_time_points();
    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points -1 ; i++) {
        m_forceofinfection +=
            parameters.get<TransmissionProbabilityOnContact>() *
            parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(m_transitions.get_last_time())(0, 0) *
            ((parameters.get<TransitionProbabilities>()[0] * 
                parameters.get<TransitionDistributions>()[1].Distribution(m_transitions.get_time(num_time_points - 1 - i)) +
              (1 - parameters.get<TransitionProbabilities>()[0]) *
                  parameters.get<TransitionDistributions>()[2].Distribution(m_transitions.get_time(num_time_points - 1  - i))) *
                 m_transitions[i+1][Eigen::Index(InfectionTransitions::ExposedToInfectedNoSymptoms)] *
                 parameters.get<RelativeTransmissionNoSymptoms>() +
             (parameters.get<TransitionProbabilities>()[1] * 
                parameters.get<TransitionDistributions>()[3].Distribution(m_transitions.get_time(num_time_points - 1 - i)) +
              (1 - parameters.get<TransitionProbabilities>()[1]) *
                  parameters.get<TransitionDistributions>()[4].Distribution(
                      m_transitions.get_time(num_time_points - 1 - i))) *
                 m_transitions[i+1][Eigen::Index(InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms)] *
                 parameters.get<RiskOfInfectionFromSymptomatic>());
    }
    printf("Force of infection tmp: %f\n", m_forceofinfection);
    m_forceofinfection = m_dt / ((ScalarType)m_N - m_SECIR.get_last_value()[Eigen::Index(InfectionState::Dead)]) * m_forceofinfection;
}

// define function that defines general scheme to compute flow
void Model::compute_flow(int idx_InfectionTransitions, ScalarType TransitionProbability, int idx_TransitionDistribution)
{
    // allowed just for idx_InfectionTransitions>=1
    ScalarType sum = 0;
    Eigen::Index calc_time_index  = (Eigen::Index)std::ceil(parameters.get<TransitionDistributions>()[idx_TransitionDistribution].get_xright() /
                                    m_dt);
    Eigen::Index num_time_points = m_transitions.get_num_time_points();
    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points-1; i++) {
        // forward difference scheme to get always tau>0
        sum += (parameters.get<TransitionDistributions>()[idx_TransitionDistribution].Distribution(
                    m_transitions.get_time(num_time_points - i)) -
                parameters.get<TransitionDistributions>()[idx_TransitionDistribution].Distribution(
                    m_transitions.get_time(num_time_points - 1 - i))) /
               m_dt * m_transitions[i+1][Eigen::Index(idx_InfectionTransitions - 1)];// i oder i+1?? siehe auch Formel -> Martin
    }
    
    m_transitions.get_last_value()[Eigen::Index(idx_InfectionTransitions)] =
        (-1) * TransitionProbability * m_dt * sum;
}

void Model::compute_totaldeaths()
{
    Eigen::Index num_time_points = m_SECIR.get_num_time_points();
    m_SECIR.get_last_value()[Eigen::Index(InfectionState::Dead)] = 
        m_SECIR[num_time_points - 2][Eigen::Index(InfectionState::Dead)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::InfectedCriticalToDead)];
}


// function that actually computes size of compartments from flows (all but S and R)
void Model::get_size_of_compartments(int idx_IncomingFlow, ScalarType TransitionProbability,
                                     int idx_TransitionDistribution, int idx_TransitionDistributionToRecov,
                                     int idx_InfectionState)
{
    ScalarType sum = 0;
    ScalarType calc_time      = parameters.get<TransitionDistributions>()[idx_TransitionDistributionToRecov].get_xright();
   
    if (parameters.get<TransitionDistributions>()[idx_TransitionDistribution].get_xright() > calc_time ) {
            calc_time = parameters.get<TransitionDistributions>()[idx_TransitionDistribution].get_xright();
    } 
    
    Eigen::Index calc_time_index         = (Eigen::Index)std::ceil(calc_time/ m_dt); //need calc_time/dt timesteps in sum
    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    for (Eigen::Index i = num_time_points -1 - calc_time_index; i < num_time_points-1; i++) {
        sum +=
            (TransitionProbability * parameters.get<TransitionDistributions>()[idx_TransitionDistribution].Distribution(
                                         m_transitions.get_time(num_time_points - 1 - i)) +
             (1 - TransitionProbability) *
                 parameters.get<TransitionDistributions>()[idx_TransitionDistributionToRecov].Distribution(
                     m_transitions.get_time(num_time_points - 1 - i))) *
            m_transitions[i+1][Eigen::Index(idx_IncomingFlow)]; // i oder i+1?? siehe auch Formel -> Martin
    }

    m_SECIR.get_last_value()[Eigen::Index(idx_InfectionState)] = m_dt * sum;
}

// do simulation
void Model::simulate(int t_max)
{
    std::cout << t_max << "\n";
    // compute S(0) and phi(0) as initial values for discretization scheme//eig müsste ich hier phi(-1) berechnen und S(0) oder? Sonst passt das ja alles gar nicht mit unserem Schema zusammen
    update_forceofinfection();
    m_SECIR[Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] =  
                    m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::SusceptibleToExposed)] /
                    m_forceofinfection;

    // for every time step:
    while ((int) m_transitions.get_last_time() < t_max) {
        std::cout << "time: " << m_transitions.get_last_time() << "\n";
        ScalarType last_time = m_transitions.get_last_time();
        m_transitions.add_time_point(last_time + m_dt);

        Eigen::Index idx = m_transitions.get_num_time_points() - 1;
        std::cout << "index time step: " << idx << "\n";
        // compute_S:
        // compute S(t_{idx}) 
        update_susceptibles();
       
        // compute flows:
        // compute sigma_E^C(t_{n+1}) and store it
        compute_flow(1, 1, 0);

        // compute sigma_C^I(t_{n+1}) and store it
        compute_flow(2, parameters.get<TransitionProbabilities>()[0], 1); //
    
        //sigma_I^H(t_{n+1})
        compute_flow(3, parameters.get<TransitionProbabilities>()[1], 3);
        //sigma_H^U(t_{n+1})
        compute_flow(4, parameters.get<TransitionProbabilities>()[2], 5);
        //sigma_U^D(t_{n+1})
        compute_flow(5, parameters.get<TransitionProbabilities>()[3], 7);

        // Compute Compartments and phi:
        // compute D
        compute_totaldeaths();

        // compute_phi(t):
        update_forceofinfection();
        
        //calculate sigma_S^E with updatet force of infection
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::SusceptibleToExposed)] =
             m_forceofinfection * m_SECIR.get_last_value()[Eigen::Index(InfectionState::Susceptible)];
        
        // compute remaining compartments from flows //TODO: Why is I(0) so groß? Fehler hier?
        // E, use dummy transitiontorecov
        get_size_of_compartments(1, parameters.get<TransitionProbabilities>()[0], 0, 0, 1);
        // C
        get_size_of_compartments(2, parameters.get<TransitionProbabilities>()[1], 1, 2, 2);
        // I
        get_size_of_compartments(3, parameters.get<TransitionProbabilities>()[2], 3, 4, 3);
        // H
        get_size_of_compartments(4, parameters.get<TransitionProbabilities>()[3], 5, 6, 4);
        // U
        get_size_of_compartments(5, parameters.get<TransitionProbabilities>()[4], 7, 8, 5);
        // R
        // TODO: needs to be updated with computation by flows
        Eigen::Index num_points_R = m_SECIR.get_num_time_points();
        m_SECIR.get_last_value()[Eigen::Index(6)] =
            m_N - m_SECIR[num_points_R - 1][Eigen::Index(0)] - m_SECIR[num_points_R - 1][Eigen::Index(1)] -
            m_SECIR.get_last_value()[Eigen::Index(2)] - m_SECIR[num_points_R - 1][Eigen::Index(3)] -
            m_SECIR.get_last_value()[Eigen::Index(4)] - m_SECIR.get_last_value()[Eigen::Index(5)] -
            m_SECIR.get_last_value()[Eigen::Index(7)];
    }       
}

void Model::print_transitions() const
{
    
    // print transitions after simulation
    std::cout << "# time  |  S -> E  |  E - > C  |  C -> I  |  I -> H  |  H -> U  |  U -> D" << std::endl;
    Eigen::Index num_points = m_transitions.get_num_time_points();
    for (Eigen::Index i = 0; i < num_points; ++i) {
        std::cout << m_transitions.get_time(i) << "      |  " << m_transitions[i][0] << "  |  "
                  << m_transitions[i][1] << "  |  " << m_transitions[i][2] << "  |  " << m_transitions[i][3]
                  << "  |  " << m_transitions[i][4] << "  |  " << m_transitions[i][5]
                  << std::endl; //[Eigen::Index(InfectionState::S)] << std::endl;
    }   
}

void Model::print_compartments() const
{
    
    // print transitions after simulation
    std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
    Eigen::Index num_points = m_SECIR.get_num_time_points();
    for (Eigen::Index i = 0; i < num_points; ++i) {
        std::cout << m_SECIR.get_time(i) << "      |  " << m_SECIR[i][0] << "  |  "
                  << m_SECIR[i][1] << "  |  " << m_SECIR[i][2] << "  |  " << m_SECIR[i][3]
                  << "  |  " << m_SECIR[i][4] << "  |  " << m_SECIR[i][5]
                  << "  |  " << m_SECIR[i][6] << "  |  " << m_SECIR[i][7] << "  |  "
                  << std::endl; //[Eigen::Index(InfectionState::S)] << std::endl;
    }   
}

} // namespace isecir
} // namespace mio
