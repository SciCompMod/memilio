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
    , m_SECIR{TimeSeries<ScalarType>(8)}
    , m_dt{dt_init}
    , m_N{N_init}
{
    m_SECIR.add_time_point(0);
    m_SECIR[Eigen::Index(0)][Eigen::Index(InfectionState::Dead)] = Dead0;
}

// define helper functions later needed for simulation
void Model::update_susceptibles()
{
    // using number of susceptibles from previous time step and force of infection from previous time step: compute current number of susceptibles
    m_susceptibles = m_susceptibles / (1 + m_dt * m_forceofinfection);

    // store susceptibles in m_SECIR
    m_SECIR[m_SECIR.get_num_time_points() - 1][Eigen::Index(InfectionState::Susceptible)] = m_susceptibles;
}

void Model::update_forceofinfection()
{
    // determine force of infection for the current last timepoint in m_transitions
    m_forceofinfection = 0;
    ScalarType calc_time      = 0;
    // determine the relevant calculation area; just in the support of one of the transition distributions

    // berechne hier calc-time, also wie viele zeitschritte wir zurückgehen müssen; hängt von xright ab
    // kann man hier nicht auch eine Funktion nehmen, die das Maximum berechnet für eine Menge von Indizes? Gibt es sowas schon?
    // TODO: when computing calc_time via this loop we obtain calc_time = 0 although it should be 1
    for (int i = 1; i < 5; i++) {
        if (parameters.get<TransitionDistributions>()[i].get_xright() < calc_time) {
            std::cout << "i" << i << "\n";
            calc_time = parameters.get<TransitionDistributions>()[i].get_xright();
        }
    }
    // set dummy calc_time here, need to compute right calc_time, see above
    calc_time = parameters.get<TransitionDistributions>()[0].get_xright();
    std::cout << "Calc_time: " << calc_time << "\n";
    //determine corresponding indices
    //Achtung: Parameter müssten zum teil noch abhängig von tau und t sein, das ist in parameters aber noch nicht
    calc_time                    = (int)std::ceil((ScalarType)calc_time / m_dt); //need calc_time timesteps in sum
    //std::cout << "Calc_time: " << calc_time << "\n";
    Eigen::Index num_time_points = m_transitions.get_num_time_points();
    for (Eigen::Index i = num_time_points - (Eigen::Index)calc_time; i < num_time_points ; i++) {
        // contactmatrix returns 0, i dont know why, use dummy to find other errors if function
        // printf("contacts: %f \n", parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(m_transitions.get_last_time())(0, 0));
        // parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(m_transitions.get_last_time())(0, 0)
        ScalarType contactmatrixdummy = 1;
    
        m_forceofinfection +=
            parameters.get<TransmissionProbabilityOnContact>() *
            contactmatrixdummy *
            ((parameters.get<TransitionProbabilities>()[0] * parameters.get<TransitionDistributions>()[1].Distribution(
                                                                 m_transitions.get_time(num_time_points - 1 - i)) +
              (1 - parameters.get<TransitionProbabilities>()[0]) *
                  parameters.get<TransitionDistributions>()[2].Distribution(
                      m_transitions.get_time(num_time_points - 1 - i))) *
                 m_transitions[i][Eigen::Index(InfectionTransitions::ExposedToInfectedNoSymptoms)] *
                 parameters.get<RelativeTransmissionNoSymptoms>() +
             (parameters.get<TransitionProbabilities>()[1] * parameters.get<TransitionDistributions>()[3].Distribution(
                                                                 m_transitions.get_time(num_time_points - 1 - i)) +
              (1 - parameters.get<TransitionProbabilities>()[1]) *
                  parameters.get<TransitionDistributions>()[4].Distribution(
                      m_transitions.get_time(num_time_points - 1 - i))) *
                 m_transitions[i][Eigen::Index(InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms)] *
                 parameters.get<RiskOfInfectionFromSymptomatic>());
    }
    printf("Force of infection tmp: %f\n", m_forceofinfection);
    //Achtung:hier muss noch D(t_n) rein
    m_forceofinfection = m_dt / ((ScalarType)m_N) * m_forceofinfection;
}

// define function that defines general scheme to compute flow
void Model::compute_flow(int idx_InfectionTransitions, ScalarType TransitionProbability, int idx_TransitionDistribution)
{
    // std::cout << "is in compute_flow \n";
    // allowed just for idx_InfectionTransitions>=1
    ScalarType sum = 0;
    int calc_time  = (int)std::ceil(parameters.get<TransitionDistributions>()[idx_TransitionDistribution].get_xright() /
                                    (double)m_dt);
    Eigen::Index num_time_points = m_transitions.get_num_time_points();
    // std::cout << "xright " << parameters.get<TransitionDistributions>()[idx_TransitionDistribution].get_xright() << "\n";
    // std::cout << "time -calc " << num_time_points - calc_time << "\n";
    // std::cout << "time " << num_time_points << "\n";
    for (Eigen::Index i = num_time_points - calc_time; i < num_time_points; i++) {
        // std::cout << "idx " << i << "\n";

        // forward difference scheme to get always tau>0
        sum += (parameters.get<TransitionDistributions>()[idx_TransitionDistribution].Distribution(
                    m_transitions.get_time(num_time_points - i - 1)) -
                parameters.get<TransitionDistributions>()[idx_TransitionDistribution].Distribution(
                    m_transitions.get_time(num_time_points - i))) /
               m_dt * m_transitions[i][Eigen::Index(idx_InfectionTransitions - 1)];
    }
    // Use m_transitions[num_time_points-1] because we are accessing via index
    m_transitions[num_time_points - 1][Eigen::Index(idx_InfectionTransitions)] =
        (-1) * TransitionProbability * m_dt * sum;
}

void Model::compute_totaldeaths()
{
    m_SECIR.add_time_point(m_SECIR.get_last_time() + m_dt);

    Eigen::Index num_time_points = m_SECIR.get_num_time_points();

    m_SECIR[num_time_points - 1][Eigen::Index(InfectionState::Dead)] =
        m_SECIR[num_time_points - 2][Eigen::Index(InfectionState::Dead)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::InfectedCriticalToDead)];

    /*std::cout << "Computed sum for deaths"
              << "\n";
    // std::cout << "m_SECIR" << m_SECIR << "\n";

    // std::cout << idx_InfectionState << "\n";
    Eigen::Index num_points = m_SECIR.get_num_time_points();
    std::cout << "num_points: " << num_points << "\n";
    for (Eigen::Index i = 0; i < num_points; ++i) {
        std::cout << "i: " << i << "\n";
        printf(" %.9f %.9f \n", m_SECIR.get_time(i), m_SECIR[i][Eigen::Index(InfectionState::Dead)]);
    }

    //m_SECIR[m_SECIR.get_num_time_points()][Eigen::Index(InfectionState::Dead)] = m_dt * sum;
    std::cout << "Updated deaths"
              << "\n";*/
}

// function that actually computes flows for remaining transitions
// TimeSeries<ScalarType> get_flows(int t_max)
// {}
// sigma_S^E, sigma_E^C and sigma_C^I are known from update_forceofinfection
// compute remaining flows here

// define function that defines general scheme to compute compartments from flows

// fucntion that actually computes size of compartments from flows (all but S and R)
void Model::get_size_of_compartments(int idx_IncomingFlow, ScalarType TransitionProbability,
                                     int idx_TransitionDistribution, int idx_TransitionDistributionToRecov,
                                     int idx_InfectionState)
{
    ScalarType sum = 0;

    // TODO: calculate calc_time for right TransitionDistributions, need to iterate over right indices
    // // berechne hier calc-time, also wie viele zeitschritte wir zurückgehen müssen; hängt von xright ab
    // for (int i = 1; i < 5; i++) {
    //     if (parameters.get<TransitionDistributions>()[i].get_xright()<calc_time) {
    //         std::cout << "i" << i << "\n";
    //         calc_time = parameters.get<TransitionDistributions>()[i].get_xright();
    // }
    // }
    // this is just a dummy calc_time, see above, tbd
    Eigen::Index calc_time                    = (Eigen::Index)std::ceil((ScalarType)parameters.get<TransitionDistributions>()[idx_TransitionDistribution].get_xright()/ m_dt); //need calc_time timesteps in sum
    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    for (Eigen::Index i = num_time_points - calc_time; i < num_time_points; i++) {
        sum +=
            (TransitionProbability * parameters.get<TransitionDistributions>()[idx_TransitionDistribution].Distribution(
                                         m_transitions.get_time(num_time_points - i)) +
             (1 - TransitionProbability) *
                 parameters.get<TransitionDistributions>()[idx_TransitionDistributionToRecov].Distribution(
                     m_transitions.get_time(num_time_points - i))) *
            m_transitions[i][Eigen::Index(idx_IncomingFlow)];
    }

    Eigen::Index num_time_points_SECIR = m_SECIR.get_num_time_points();

    // for (Eigen::Index i = 0; i < num_time_points_SECIR; ++i) {
    //     printf(" %.9f %.9f \n", m_SECIR.get_time(i), m_SECIR[i][Eigen::Index(idx_InfectionState)]);
    // }
    m_SECIR[num_time_points_SECIR - 1][Eigen::Index(idx_InfectionState)] = m_dt * sum;
}

// do simulation
void Model::simulate(int t_max)
{
    std::cout << t_max << "\n";
    // compute S(0) and phi(0) as initial values for discretization scheme
    update_forceofinfection();
    m_susceptibles = m_transitions[m_transitions.get_num_time_points() - 1]
                                  [Eigen::Index(InfectionTransitions::SusceptibleToExposed)] /
                     m_forceofinfection;

    // for every time step:
    while (m_transitions.get_last_time() < t_max) {
        std::cout << "time: " << m_transitions.get_last_time() << "\n";
        ScalarType last_time = m_transitions.get_last_time();
        m_transitions.add_time_point(last_time + m_dt);

        Eigen::Index idx = m_transitions.get_num_time_points() - 1;
        std::cout << "index: " << idx << "\n";
        // compute_S:
        // compute S(t_{idx}) and store it in m_susceptibles
        update_susceptibles();
        std::cout << "Updated susceptibles"
                  << "\n";
        std::printf("Susceptibles: %f \n", m_susceptibles);
        // compute_phi:
        // compute sigma_E^C(t_{n+1}) and store it
        compute_flow(1, 1, 0);
        // std::cout << "Computed flow 1"
        //           << "\n";
        // compute sigma_C^I(t_{n+1}) and store it
        compute_flow(2, parameters.get<TransitionProbabilities>()[0], 1); //
        // std::cout << "Computed flow 2"
        //           << "\n";
        // with this: compute phi(t_{n+1})
        update_forceofinfection();
        std::cout << "Updated forceofinfection"
                  << "\n";
        std::printf("Force of infection: %f \n", m_forceofinfection);
        //calculate sigma_S^E
        m_transitions[idx][Eigen::Index(InfectionTransitions::SusceptibleToExposed)] =
            m_forceofinfection * m_susceptibles;

        // compute remaining flows
        //sigma_I^H(t_{n+1})
        compute_flow(3, parameters.get<TransitionProbabilities>()[1], 3);
        //sigma_H^U(t_{n+1})
        compute_flow(4, parameters.get<TransitionProbabilities>()[2], 5);
        //sigma_U^D(t_{n+1})
        compute_flow(5, parameters.get<TransitionProbabilities>()[3], 7);
        std::cout << "Computed all flows"
                  << "\n";
        // compute D
        compute_totaldeaths();
        std::cout << "Computed deaths"
                  << "\n";
    }

    // at the end of simulation:
    // compute remaining compartments from flows
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

    std::cout << "almost done ... \n";
    // R
    // needs to be updated with computation by flows
    Eigen::Index num_points_R = m_SECIR.get_num_time_points();
    m_SECIR[num_points_R - 1][Eigen::Index(6)] =
        m_N - m_SECIR[num_points_R - 1][Eigen::Index(0)] - m_SECIR[num_points_R - 1][Eigen::Index(1)] -
        m_SECIR[num_points_R - 1][Eigen::Index(2)] - m_SECIR[num_points_R - 1][Eigen::Index(3)] -
        m_SECIR[num_points_R - 1][Eigen::Index(4)] - m_SECIR[num_points_R - 1][Eigen::Index(5)] -
        m_SECIR[num_points_R - 1][Eigen::Index(7)];
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
