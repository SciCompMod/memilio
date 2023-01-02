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

void Model::update_susceptibles()
{
    Eigen::Index num_time_points = m_SECIR.get_num_time_points();
    // using number of susceptibles from previous time step and force of infection from previous time step:
    // compute current number of susceptibles and store susceptibles in m_SECIR
    m_SECIR.get_last_value()[Eigen::Index(InfectionState::Susceptible)] =
        m_SECIR[num_time_points - 2][Eigen::Index(InfectionState::Susceptible)] / (1 + m_dt * m_forceofinfection);
}

void Model::update_forceofinfection()
{
    // determine force of infection for the current last point of time in m_transitions,
    m_forceofinfection   = 0;
    ScalarType calc_time = 0;

    // determine the relevant calculation area = union of the supports of the relevant transition distributions
    // TODO: evtl gibt es eine Funktion, die das maximum aus verschiedenen Werten nimmt?
    for (int i = 1; i < 5; i++) {
        if (parameters.get<TransitionDistributions>()[i].get_xright() > calc_time) {
            calc_time = parameters.get<TransitionDistributions>()[i].get_xright();
        }
    }

    //corresponding index
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / m_dt); //need calc_time timesteps in sum

    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    // TODO: Parameter müssten zum teil noch abhängig von tau und t sein, das ist in parameters aber noch nicht
    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {
        m_forceofinfection +=
            parameters.get<TransmissionProbabilityOnContact>() *
            parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(m_transitions.get_last_time())(0, 0) *
            ((parameters.get<TransitionProbabilities>()[0] *
                  parameters.get<TransitionDistributions>()[1].Distribution((num_time_points - 1 - i) * m_dt) +
              (1 - parameters.get<TransitionProbabilities>()[0]) *
                  parameters.get<TransitionDistributions>()[2].Distribution((num_time_points - 1 - i) * m_dt)) *
                 m_transitions[i + 1][Eigen::Index(InfectionTransitions::ExposedToInfectedNoSymptoms)] *
                 parameters.get<RelativeTransmissionNoSymptoms>() +
             (parameters.get<TransitionProbabilities>()[1] *
                  parameters.get<TransitionDistributions>()[3].Distribution((num_time_points - 1 - i) * m_dt) +
              (1 - parameters.get<TransitionProbabilities>()[1]) *
                  parameters.get<TransitionDistributions>()[4].Distribution((num_time_points - 1 - i) * m_dt)) *
                 m_transitions[i + 1][Eigen::Index(InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms)] *
                 parameters.get<RiskOfInfectionFromSymptomatic>());

        m_forceofinfection = m_dt / ((ScalarType)m_N - m_SECIR.get_last_value()[Eigen::Index(InfectionState::Dead)]) *
                             m_forceofinfection;
    }
}

// define function that defines general scheme to compute flow
void Model::compute_flow(int idx_InfectionTransitions, ScalarType TransitionProbability, int idx_TransitionDistribution)
{ // allowed just for idx_InfectionTransitions>=1

    ScalarType sum = 0;
    // We have Transitiondistribution(m_dt*i)=0 for all i>= calc_time_index
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(
        parameters.get<TransitionDistributions>()[idx_TransitionDistribution].get_xright() / m_dt);
    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {
        // (num_time_points - 1 - i)* m_dt is the time, the individuals already been infected
        ScalarType infection_age = (num_time_points - 1 - i) * m_dt;

        // backward difference scheme to approximate first derivative
        sum +=
            (parameters.get<TransitionDistributions>()[idx_TransitionDistribution].Distribution(infection_age) -
             parameters.get<TransitionDistributions>()[idx_TransitionDistribution].Distribution(infection_age - m_dt)) /
            m_dt * m_transitions[i + 1][Eigen::Index(idx_InfectionTransitions - 1)];
        // TODO: müssen überlegen ob i oder i+1 für Zeitindex in m_transitions?? siehe auch Formel -> Martin
    }

    m_transitions.get_last_value()[Eigen::Index(idx_InfectionTransitions)] = (-1) * TransitionProbability * m_dt * sum;
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

    // deermine relevant calculation area and corresponding index
    ScalarType calc_time = parameters.get<TransitionDistributions>()[idx_TransitionDistributionToRecov].get_xright();
    if (parameters.get<TransitionDistributions>()[idx_TransitionDistribution].get_xright() > calc_time) {
        calc_time = parameters.get<TransitionDistributions>()[idx_TransitionDistribution].get_xright();
    }
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / m_dt);

    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {
        sum +=
            (TransitionProbability * parameters.get<TransitionDistributions>()[idx_TransitionDistribution].Distribution(
                                         (num_time_points - 1 - i) * m_dt) +
             (1 - TransitionProbability) *
                 parameters.get<TransitionDistributions>()[idx_TransitionDistributionToRecov].Distribution(
                     (num_time_points - 1 - i) * m_dt)) *
            m_transitions[i + 1][Eigen::Index(idx_IncomingFlow)];
        // TODO: i oder i+1?? siehe auch Formel -> Martin
    }

    m_SECIR.get_last_value()[Eigen::Index(idx_InfectionState)] = m_dt * sum;
}

// do simulation
void Model::simulate(int t_max)
{
    std::cout << "Starting simulation:  \n";
    // compute S(0) and phi(0) as initial values for discretization scheme
    update_forceofinfection();
    m_SECIR[Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] =
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::SusceptibleToExposed)] / m_forceofinfection;

    //calculate compartment sizes for t=0
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
    m_SECIR.get_last_value()[Eigen::Index(6)] =
        m_N - m_SECIR.get_last_value()[Eigen::Index(0)] - m_SECIR.get_last_value()[Eigen::Index(1)] -
        m_SECIR.get_last_value()[Eigen::Index(2)] - m_SECIR.get_last_value()[Eigen::Index(3)] -
        m_SECIR.get_last_value()[Eigen::Index(4)] - m_SECIR.get_last_value()[Eigen::Index(5)] -
        m_SECIR.get_last_value()[Eigen::Index(7)];

    // for every time step:
    while ((int)m_transitions.get_last_time() < t_max) {

        m_transitions.add_time_point(m_transitions.get_last_time() + m_dt);
        std::cout << "time: " << m_transitions.get_last_time() << "\n";
        m_SECIR.add_time_point(m_SECIR.get_last_time() + m_dt);

        // compute_S:
        update_susceptibles();

        //calculate sigma_S^E with force of infection from previous time step und susceptibles from current time step
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::SusceptibleToExposed)] =
            m_forceofinfection * m_SECIR.get_last_value()[Eigen::Index(InfectionState::Susceptible)];

        // compute flows:
        //TODO: Hier wollten wir auch die Flows nach R berechnen
        // compute sigma_E^C
        compute_flow(1, 1, 0);
        // compute sigma_C^I
        compute_flow(2, parameters.get<TransitionProbabilities>()[0], 1);
        //sigma_I^H
        compute_flow(3, parameters.get<TransitionProbabilities>()[1], 3);
        //sigma_H^U
        compute_flow(4, parameters.get<TransitionProbabilities>()[2], 5);
        //sigma_U^D
        compute_flow(5, parameters.get<TransitionProbabilities>()[3], 7);

        // compute D
        compute_totaldeaths();

        // compute phi(t) (only used for calculation of S in the next timestep!):
        update_forceofinfection();

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
        // R
        // TODO: needs to be updated with computation by flows
        m_SECIR.get_last_value()[Eigen::Index(6)] =
            m_N - m_SECIR.get_last_value()[Eigen::Index(0)] - m_SECIR.get_last_value()[Eigen::Index(1)] -
            m_SECIR.get_last_value()[Eigen::Index(2)] - m_SECIR.get_last_value()[Eigen::Index(3)] -
            m_SECIR.get_last_value()[Eigen::Index(4)] - m_SECIR.get_last_value()[Eigen::Index(5)] -
            m_SECIR.get_last_value()[Eigen::Index(7)];
    }
    std::cout << "Finished simulation successfully.  \n";
}

/*TODO: maybe add a print function to Timeseries class instead of writing for each timeseries a new one
Als Übergabeparameter könte man auch String Ueberschrift oder sowas machen, damit die Form gut lesbar bleibt.*/
void Model::print_transitions() const
{

    // print transitions after simulation
    std::cout << "# time  |  S -> E  |  E - > C  |  C -> I  |  I -> H  |  H -> U  |  U -> D  " << std::endl;
    Eigen::Index num_points = m_transitions.get_num_time_points();
    for (Eigen::Index i = 0; i < num_points; ++i) {
        std::cout << m_transitions.get_time(i) << "      |  " << m_transitions[i][0] << "  |  " << m_transitions[i][1]
                  << "  |  " << m_transitions[i][2] << "  |  " << m_transitions[i][3] << "  |  " << m_transitions[i][4]
                  << "  |  " << m_transitions[i][5] << "  |  " << std::endl;
    }
}

void Model::print_compartments() const
{

    // print transitions after simulation
    std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
    Eigen::Index num_points = m_SECIR.get_num_time_points();
    for (Eigen::Index i = 0; i < num_points; ++i) {
        std::cout << m_SECIR.get_time(i) << "      |  " << m_SECIR[i][0] << "  |  " << m_SECIR[i][1] << "  |  "
                  << m_SECIR[i][2] << "  |  " << m_SECIR[i][3] << "  |  " << m_SECIR[i][4] << "  |  " << m_SECIR[i][5]
                  << "  |  " << m_SECIR[i][6] << "  |  " << m_SECIR[i][7] << "  |  " << std::endl;
    }
}

} // namespace isecir
} // namespace mio
