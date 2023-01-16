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

#ifndef IDESECIR_MODEL_H
#define IDESECIR_MODEL_H

#include "ide_secir/parameters.h"
#include "ide_secir/infection_state.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isecir
{
class Model
{
    /* TODO: 
    - in parameters muessen einige Parameter noch von tau abhaengig gemacht werden (ContactMatrix)
    - die Dokumentationen der Funktionen sind unvollstaendig
    - i oder i+1 in numerischer Integration
    -ueberlegen ob eine allgemeine print -funktion fuer timeSeries Sinn ergibt. aktuell haben wir hier 2 mal und in SEIR nochmal dieselbe print-Funktion
    - wir sollten eine "constraint check"- Funktion in Parameters schreiben, die zB prüft od Probability C->I =1- Prob C->R ist.
    In der könnte man auch S->E auf 1 setzen oder die dummys generell auf 0/1/NaN. 
    evtl ist dafür dieses Array in Infection state gar nicht so schlecht, weil man darüber sehen kann, wie viele flows von einem Kompartiment ausgehen. die summe der zugehoerigen Wahrscheinlichkeiten muss 1 sein
    - Evtl. allgemeinen Integrator als non-standrd diff. scheme anstatt aktueller Approximation
    */
    using Pa = ParametersBase;

public:
    /**
    * @brief Constructor to create an IDE SECIR model.
    *
    * @param[in, out] init TimeSeries with the initial values of the number of individuals,
    *   which transition according to the transitions defined in infection_state.h at associated initial times.
    *   The time steps should be equidistant and equal to the time step used for the simulation. 
    *   A certain history of time steps and values for the transitions is needed. 
    *   A warning is displayed if the condition is violated.
    *   The last time point in this vector should be a time 0.
    * @param[in] dt_init The size of the time step used for numerical simulation.
    * @param[in] N_init The population of the considered region. 
    * @param[in, out] Parameterset_init used Parameters for simulation. 
    */
    Model(TimeSeries<ScalarType>&& init, ScalarType dt_init, size_t N_init, size_t Dead0, Pa Parameterset_init = Pa());

    /**
    * @brief Simulate the evolution of infection numbers with the given IDE SECIR model.
    *
    * The simulation is performed by solving the underlying model equation numerically. 
    * Here, an integro-differential equation is to be solved. The model parameters and the initial data are used.
    *
    * @param[in] t_max Last simulation day. 
    *   If the last point of time of the initial TimeSeries was 0, the simulation will be executed for t_max days.
    * @return The result of the simulation, stored in a TimeSeries with simulation time and 
    *       associated number of susceptibles.
    */
    void simulate(int t_max);
    // Used Parameters for the simulation.
    Pa parameters{};

    /**
     * @brief print the transition part of the simulation result.
     * 
     * The TimeSeries m_transitions with initial values used for the simulation and calculated transitions by the 
     * simulation are printed. 
     */
    void print_transitions() const;

    /**
     * @brief print the simulated numbers of individuals in each compartment for each time step.
     * 
     * The TimeSeries m_SECIR with simulated numbers of individuals in each compartment for each time step are printed. 
     */
    void print_compartments() const;

private:
    void update_susceptibles();
    void update_forceofinfection();
    void compute_flow(int idx_InfectionTransitions, Eigen::Index idx_IncomingFlow);
    void update_flows();
    void compute_totaldeaths();
    void compute_recovered();

    /**
     * @brief Get the size of the compartment specified in idx_InfectionState at the current last time in m_SECIR.
     * 
     * @param[in] idx_InfectionState specifies the considered InfectionState
     * @param[in] idx_IncomingFlow specifies the index of the infoming flow to InfectionState in m_transitions. 
     * @param[in] idx_TransitionDistribution1 specifies the index of the first relevant transitiondistribution, 
     *              related to a flow from the considered InfectionState to any other State.
     *              This index is also used for related Probability.
     * @param[in] idx_TransitionDistribution2 specifies the index of the second relevant transitiondistribution, 
     *              related to a flow from the considered InfectionState to any other State (in most cases to Recovered). 
     * @param[in] transitionprobability2 transitionsprobability related to idx_TransitionDistribution2.
     *              This is just an extra parameter to give the opportunity that a compartment only have one outgoing 
     *              flow (eg as E). In this case this probability has to be set to 0.
     *              If the InfectionState has two outgoing flows, one could also get this related probability 
     *              via idx_TransitionDistribution2.
     */
    void get_size_of_compartments(Eigen::Index idx_InfectionState, Eigen::Index idx_IncomingFlow,
                                  int idx_TransitionDistribution1, int idx_TransitionDistribution2);
    void update_compartments_ECIHU();

    // TimeSeries containing points of time and the corresponding number of transitions.
    TimeSeries<ScalarType> m_transitions;
    // TimeSeries containing points of time and the corresponding number of people in defined Infectionstates.
    TimeSeries<ScalarType> m_SECIR;
    /*Attention: m_SECIR and m_transitions do not necessarily have the same number of time points*/

    // Force of infection term needed for numerical scheme, corresponds to phi
    ScalarType m_forceofinfection{0};

    // Timestep used for simulation.
    ScalarType m_dt{0};
    // Population of the considered region.
    size_t m_N{0};
};

} // namespace isecir
} // namespace mio

#endif // IDESECIR_MODEL_H
