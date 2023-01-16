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
#include "memilio/utils/logging.h"
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
    if (!((int)m_transitions.get_num_elements() == (int)InfectionTransitions::Count)) {
        log_error("Initialization failed. Number of elements in init does not match the required number.");
    }
    m_SECIR.add_time_point(0);
    m_SECIR[Eigen::Index(0)][Eigen::Index(InfectionState::Dead)] = Dead0;
}

void Model::compute_susceptibles()
{
    Eigen::Index num_time_points = m_SECIR.get_num_time_points();
    // using number of susceptibles from previous time step and force of infection from previous time step:
    // compute current number of susceptibles and store susceptibles in m_SECIR
    m_SECIR.get_last_value()[Eigen::Index(InfectionState::Susceptible)] =
        m_SECIR[num_time_points - 2][Eigen::Index(InfectionState::Susceptible)] / (1 + m_dt * m_forceofinfection);
}

void Model::update_forceofinfection()
{
    m_forceofinfection = 0;

    // determine the relevant calculation area = union of the supports of the relevant transition distributions
    ScalarType calc_time = std::max(
        {parameters.get<TransitionDistributions>()[(int)InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms]
             .get_xright(),
         parameters.get<TransitionDistributions>()[(int)InfectionTransitions::InfectedNoSymptomsToRecovered]
             .get_xright(),
         parameters.get<TransitionDistributions>()[(int)InfectionTransitions::InfectedSymptomsToInfectedSevere]
             .get_xright(),
         parameters.get<TransitionDistributions>()[(int)InfectionTransitions::InfectedSymptomsToRecovered]
             .get_xright()});

    //corresponding index
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / m_dt); //need calc_time_index timesteps in sum

    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    // TODO: Parameter müssten zum teil noch abhängig von tau und t sein, das ist in parameters aber noch nicht
    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {
        m_forceofinfection +=
            parameters.get<TransmissionProbabilityOnContact>() *
            parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(m_transitions.get_last_time())(0, 0) *
            ((parameters
                      .get<TransitionProbabilities>()[(int)InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms] *
                  parameters
                      .get<TransitionDistributions>()[(int)InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms]
                      .Distribution((num_time_points - 1 - i) * m_dt) +
              parameters.get<TransitionProbabilities>()[(int)InfectionTransitions::InfectedNoSymptomsToRecovered] *
                  parameters.get<TransitionDistributions>()[(int)InfectionTransitions::InfectedNoSymptomsToRecovered]
                      .Distribution((num_time_points - 1 - i) * m_dt)) *
                 m_transitions[i + 1][Eigen::Index(InfectionTransitions::ExposedToInfectedNoSymptoms)] *
                 parameters.get<RelativeTransmissionNoSymptoms>() +
             (parameters.get<TransitionProbabilities>()[(int)InfectionTransitions::InfectedSymptomsToInfectedSevere] *
                  parameters.get<TransitionDistributions>()[(int)InfectionTransitions::InfectedSymptomsToInfectedSevere]
                      .Distribution((num_time_points - 1 - i) * m_dt) +
              parameters.get<TransitionProbabilities>()[(int)InfectionTransitions::InfectedSymptomsToRecovered] *
                  parameters.get<TransitionDistributions>()[(int)InfectionTransitions::InfectedSymptomsToRecovered]
                      .Distribution((num_time_points - 1 - i) * m_dt)) *
                 m_transitions[i + 1][Eigen::Index(InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms)] *
                 parameters.get<RiskOfInfectionFromSymptomatic>());

        m_forceofinfection = m_dt / ((ScalarType)m_N - m_SECIR.get_last_value()[Eigen::Index(InfectionState::Dead)]) *
                             m_forceofinfection;
    }
}

void Model::compute_flow(int idx_InfectionTransitions, Eigen::Index idx_IncomingFlow)
{
    ScalarType sum = 0;

    // We have Transitiondistribution(m_dt*i)=0 for all i>= calc_time_index
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(
        parameters.get<TransitionDistributions>()[idx_InfectionTransitions].get_xright() / m_dt);
    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {
        // (num_time_points - 1 - i)* m_dt is the time, the individuals already been infected
        ScalarType infection_age = (num_time_points - 1 - i) * m_dt;

        // backward difference scheme to approximate first derivative
        sum +=
            (parameters.get<TransitionDistributions>()[idx_InfectionTransitions].Distribution(infection_age) -
             parameters.get<TransitionDistributions>()[idx_InfectionTransitions].Distribution(infection_age - m_dt)) /
            m_dt * m_transitions[i + 1][idx_IncomingFlow];
        // TODO: müssen überlegen ob i oder i+1 für Zeitindex in m_transitions?? siehe auch Formel -> Martin
    }

    m_transitions.get_last_value()[Eigen::Index(idx_InfectionTransitions)] =
        (-1) * parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * m_dt * sum;
}

void Model::flows_current_timestep()
{
    // calculate sigma_S^E with force of infection from previous time step und susceptibles from current time step
    m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::SusceptibleToExposed)] =
        m_forceofinfection * m_SECIR.get_last_value()[Eigen::Index(InfectionState::Susceptible)];
    // calculate all other flows with compute_flow
    // sigma_E^C
    compute_flow((int)InfectionTransitions::ExposedToInfectedNoSymptoms,
                 Eigen::Index(InfectionTransitions::SusceptibleToExposed));
    // sigma_C^I
    compute_flow((int)InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms,
                 Eigen::Index(InfectionTransitions::ExposedToInfectedNoSymptoms));
    // sigma_C^R
    compute_flow((int)InfectionTransitions::InfectedNoSymptomsToRecovered,
                 Eigen::Index(InfectionTransitions::ExposedToInfectedNoSymptoms));
    // sigma_I^H
    compute_flow((int)InfectionTransitions::InfectedSymptomsToInfectedSevere,
                 Eigen::Index(InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms));
    // sigma_I^R
    compute_flow((int)InfectionTransitions::InfectedSymptomsToRecovered,
                 Eigen::Index(InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms));
    // sigma_H^U
    compute_flow((int)InfectionTransitions::InfectedSevereToInfectedCritical,
                 Eigen::Index(InfectionTransitions::InfectedSymptomsToInfectedSevere));
    // sigma_H^R
    compute_flow((int)InfectionTransitions::InfectedSevereToRecovered,
                 Eigen::Index(InfectionTransitions::InfectedSymptomsToInfectedSevere));
    // sigma_U^D
    compute_flow((int)InfectionTransitions::InfectedCriticalToDead,
                 Eigen::Index(InfectionTransitions::InfectedSevereToInfectedCritical));
    // sigma_U^R
    compute_flow((int)InfectionTransitions::InfectedCriticalToRecovered,
                 Eigen::Index(InfectionTransitions::InfectedSevereToInfectedCritical));
}

void Model::compute_totaldeaths()
{
    Eigen::Index num_time_points = m_SECIR.get_num_time_points();

    m_SECIR.get_last_value()[Eigen::Index(InfectionState::Dead)] =
        m_SECIR[num_time_points - 2][Eigen::Index(InfectionState::Dead)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::InfectedCriticalToDead)];
}

void Model::compute_recovered()
{
    Eigen::Index num_time_points = m_SECIR.get_num_time_points();

    m_SECIR.get_last_value()[Eigen::Index(InfectionState::Recovered)] =
        m_SECIR[num_time_points - 2][Eigen::Index(InfectionState::Recovered)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::InfectedNoSymptomsToRecovered)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::InfectedSymptomsToRecovered)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::InfectedSevereToRecovered)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::InfectedCriticalToRecovered)];
}

void Model::get_size_of_compartments(Eigen::Index idx_InfectionState, Eigen::Index idx_IncomingFlow,
                                     int idx_TransitionDistribution1, int idx_TransitionDistribution2)
{
    ScalarType sum = 0;

    // determine relevant calculation area and corresponding index
    ScalarType calc_time =
        std::max(parameters.get<TransitionDistributions>()[idx_TransitionDistribution1].get_xright(),
                 parameters.get<TransitionDistributions>()[idx_TransitionDistribution2].get_xright());

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / m_dt);

    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {
        sum += (parameters.get<TransitionProbabilities>()[idx_TransitionDistribution1] *
                    parameters.get<TransitionDistributions>()[idx_TransitionDistribution1].Distribution(
                        (num_time_points - 1 - i) * m_dt) +
                (1 - parameters.get<TransitionProbabilities>()[idx_TransitionDistribution1]) *
                    parameters.get<TransitionDistributions>()[idx_TransitionDistribution2].Distribution(
                        (num_time_points - 1 - i) * m_dt)) *
               m_transitions[i + 1][idx_IncomingFlow];
        // TODO: i oder i+1?? siehe auch Formel -> Martin
    }

    m_SECIR.get_last_value()[idx_InfectionState] = m_dt * sum;
}

void Model::compartments_current_timestep_ECIHU()
{
    // E
    get_size_of_compartments(Eigen::Index(InfectionState::Exposed),
                             Eigen::Index(InfectionTransitions::SusceptibleToExposed),
                             (int)InfectionTransitions::ExposedToInfectedNoSymptoms, 0);
    // C
    get_size_of_compartments(Eigen::Index(InfectionState::InfectedNoSymptoms),
                             Eigen::Index(InfectionTransitions::ExposedToInfectedNoSymptoms),
                             (int)InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms,
                             (int)InfectionTransitions::InfectedNoSymptomsToRecovered);
    // I
    get_size_of_compartments(Eigen::Index(InfectionState::InfectedSymptoms),
                             Eigen::Index(InfectionTransitions::InfectedNoSymptomsToInfectedSymptoms),
                             (int)InfectionTransitions::InfectedSymptomsToInfectedSevere,
                             (int)InfectionTransitions::InfectedSymptomsToRecovered);
    // H
    get_size_of_compartments(Eigen::Index(InfectionState::InfectedSevere),
                             Eigen::Index(InfectionTransitions::InfectedSymptomsToInfectedSevere),
                             (int)InfectionTransitions::InfectedSevereToInfectedCritical,
                             (int)InfectionTransitions::InfectedSevereToRecovered);
    // U
    get_size_of_compartments(Eigen::Index(InfectionState::InfectedCritical),
                             Eigen::Index(InfectionTransitions::InfectedSevereToInfectedCritical),
                             (int)InfectionTransitions::InfectedCriticalToDead,
                             (int)InfectionTransitions::InfectedCriticalToRecovered);
}

void Model::simulate(int t_max)
{
    std::cout << "Starting simulation:  \n";
    // compute S(0) and phi(0) as initial values for discretization scheme
    update_forceofinfection();
    m_SECIR[Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] =
        m_transitions.get_last_value()[Eigen::Index(InfectionTransitions::SusceptibleToExposed)] / m_forceofinfection;

    //calculate compartment sizes for t=0
    compartments_current_timestep_ECIHU();

    //R; need an initial value for R, therefore do not calculate via compute_recovered()
    m_SECIR.get_last_value()[Eigen::Index(InfectionState::Recovered)] =
        m_N - m_SECIR.get_last_value()[Eigen::Index(InfectionState::Susceptible)] -
        m_SECIR.get_last_value()[Eigen::Index(InfectionState::Exposed)] -
        m_SECIR.get_last_value()[Eigen::Index(InfectionState::InfectedNoSymptoms)] -
        m_SECIR.get_last_value()[Eigen::Index(InfectionState::InfectedSymptoms)] -
        m_SECIR.get_last_value()[Eigen::Index(InfectionState::InfectedSevere)] -
        m_SECIR.get_last_value()[Eigen::Index(InfectionState::InfectedCritical)] -
        m_SECIR.get_last_value()[Eigen::Index(InfectionState::Dead)];

    // for every time step:
    while ((int)m_transitions.get_last_time() < t_max) {

        m_transitions.add_time_point(m_transitions.get_last_time() + m_dt);
        std::cout << "time: " << m_transitions.get_last_time() << "\n";
        m_SECIR.add_time_point(m_SECIR.get_last_time() + m_dt);

        // compute_S:
        compute_susceptibles();

        // compute flows:
        flows_current_timestep();

        // compute D
        compute_totaldeaths();

        // compute phi(t) (only used for calculation of S in the next timestep!):
        update_forceofinfection();

        // compute remaining compartments from flows
        compartments_current_timestep_ECIHU();
        compute_recovered();
    }
    std::cout << "Finished simulation successfully.  \n";
}

/*TODO: maybe add a print function to Timeseries class instead of writing for each timeseries a new one
Als Übergabeparameter könte man auch String Ueberschrift oder sowas machen, damit die Form gut lesbar bleibt.*/
void Model::print_transitions() const
{
    // print transitions after simulation
    std::cout << "# time  |  S -> E  |  E - > C  |  C -> I  |  C -> R  |  I -> H  |  I -> R  |  H -> U  |  H -> R  |  "
                 "U -> D  |  U -> R  "
              << std::endl;
    for (Eigen::Index i = 0; i < m_transitions.get_num_time_points(); ++i) {
        std::cout << m_transitions.get_time(i);
        for (Eigen::Index j = 0; j < m_transitions.get_num_elements(); ++j) {
            std::cout << "  |  " << m_transitions[i][j];
        }
        std::cout << "\n" << std::endl;
    }
}

void Model::print_compartments() const
{
    // print compartments after simulation
    std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
    for (Eigen::Index i = 0; i < m_SECIR.get_num_time_points(); ++i) {
        std::cout << m_SECIR.get_time(i);
        for (Eigen::Index j = 0; j < m_SECIR.get_num_elements(); ++j) {
            std::cout << "  |  " << m_SECIR[i][j];
        }
        std::cout << "\n" << std::endl;
    }
}

} // namespace isecir
} // namespace mio
