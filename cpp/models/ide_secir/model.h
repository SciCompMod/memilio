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
#include "memilio/config.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isecir
{
class Model
{
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
    * @param[in] Dead0 The total number of deaths at initial time 0.
    * @param[in, out] Parameterset_init used Parameters for simulation. 
    */
    Model(TimeSeries<ScalarType>&& init, ScalarType dt_init, ScalarType N_init, ScalarType Dead0,
          const Pa& Parameterset_init = Pa());

    void initialize();

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
    TimeSeries<ScalarType> const& simulate(ScalarType t_max);

    // Used Parameters for the simulation.
    Pa parameters{};

    /**
     * @brief Print the transition part of the simulation result.
     * 
     * The TimeSeries m_transitions with initial values used for the simulation and calculated transitions by the 
     * simulation are printed. 
     */
    void print_transitions() const;

    /**
     * @brief Print the simulated numbers of individuals in each compartment for each time step.
     * 
     * The TimeSeries m_SECIR with simulated numbers of individuals in each compartment for each time step are printed. 
     */
    void print_compartments() const;

    TimeSeries<ScalarType> const& get_flows();

    ScalarType get_timestep();

private:
    /**
    * @brief Computes number of Susceptibles for the current last time in m_SECIR.
    *
    * Number is computed using previous number of Susceptibles and the force of infection (also from previous timestep).
    * Number is stored at the matching index in m_SECIR.
    */
    void compute_susceptibles();

    /**
     * @brief Computes force of infection for the current last time in m_transitions.
     * 
     * Computed value is stored in m_forceofinfection.
     */
    void update_forceofinfection();

    /**
     * @brief Computes size of a flow.
     * 
     * Computes size of one flow from InfectionTransitions, specified in idx_InfectionTransitions, for the current 
     * last timevalue in m_transitions.
     *
     * @param[in] idx_InfectionTransitions Specifies the considered flow from InfectionTransitions.
     * @param[in] idx_IncomingFlow Index of the flow in InfectionTransitions, which goes to the considered starting
     *      compartment of the flow specified in idx_InfectionTransitions. Size of considered flow is calculated via 
     *      the value of this incoming flow.
     */
    void compute_flow(int idx_InfectionTransitions, Eigen::Index idx_IncomingFlow);

    /**
     * @brief Sets all required flows for the current last timestep in m_transitions.
     *
     * New values are stored in m_transitions. Most values are computed via the function compute_flow().
     * 
     */
    void flows_current_timestep();

    /**
     * @brief Computes total number of Deaths for the current last time in m_SECIR.
     * 
     * Number is stored in m_SECIR.
     *
     */
    void compute_totaldeaths();

    /**
     * @brief Computes total number of Recovered for the current last time in m_SECIR.
     * 
     * Number is stored in m_SECIR.
     *
     */
    void compute_recovered();

    /**
     * @brief Get the size of the compartment specified in idx_InfectionState at the current last time in m_SECIR.
     * 
     * Calculation is reasonable for all compartments except S, R, D. 
     * Therefore, we have alternative funtions for those compartments.
     *
     * @param[in] idx_InfectionState Specifies the considered InfectionState
     * @param[in] idx_IncomingFlow Specifies the index of the infoming flow to InfectionState in m_transitions. 
     * @param[in] idx_TransitionDistribution1 Specifies the index of the first relevant transitiondistribution, 
     *              related to a flow from the considered InfectionState to any other State.
     *              This index is also used for related Probability.
     * @param[in] idx_TransitionDistribution2 Specifies the index of the second relevant transitiondistribution, 
     *              related to a flow from the considered InfectionState to any other State (in most cases to Recovered). 
     *              Necessary related probability is calculated via 1-probability[idx_TransitionDistribution1].
     */
    void compute_compartment(Eigen::Index idx_InfectionState, Eigen::Index idx_IncomingFlow,
                             int idx_TransitionDistribution1, int idx_TransitionDistribution2);

    /**
     * @brief Sets all values of remaining compartments ECIHU for the current last timestep in m_SECIR.
     *
     * New values are stored in m_SECIR. Most values are computed via the function get_size_of_compartments().
     * 
     */
    void compartments_current_timestep_ECIHU();

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
    ScalarType m_N{0};
};

} // namespace isecir
} // namespace mio

#endif // IDESECIR_MODEL_H
