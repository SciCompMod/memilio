/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler, Lena Ploetzke, Martin J. Kuehn
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
    using ParameterSet = Parameters;

public:
    /**
    * @brief Constructor to create an IDE SECIR model.
    *
    * @param[in, out] init TimeSeries with the initial values of the number of individuals, 
    *   which transit within one timestep dt_init from one compartment to another.
    *   Possible transitions are specified in as #InfectionTransition%s.
    *   Considered points of times should have the distance dt_init and the last time point should be 0. 
    *   The time history must reach a certain point in the past so that the simulation can be performed.
    *   A warning is displayed if the condition is violated.
    * @param[in] dt_init The size of the time step used for numerical simulation.
    * @param[in] N_init The population of the considered region.
    * @param[in] Dead_before The total number of deaths at the time point - dt_init.
    * @param[in, out] Parameterset_init Used Parameters for simulation. 
    */
    Model(TimeSeries<ScalarType>&& init, ScalarType N_init, ScalarType Dead_before,
          const ParameterSet& Parameterset_init = ParameterSet());

    /**
    * @brief Checks constraints on model parameters.
    */
    void check_constraints(ScalarType dt) const
    {
        if (!(m_populations.get_num_time_points() > 0)) {
            log_error("Model construction failed. No initial time point for populations.");
        }

        for (int i = 0; i < (int)InfectionState::Count; i++) {
            if (m_populations[0][i] < 0) {
                log_error("Initialization failed. Initial values for populations are less than zero.");
            }
        }

        if (!((int)m_transitions.get_num_elements() == (int)InfectionTransition::Count)) {
            log_error(
                "Initialization failed. Number of elements in transition vector does not match the required number.");
        }

        ScalarType support_max = std::max(
            {parameters.get<TransitionDistributions>()[(int)InfectionTransition::ExposedToInfectedNoSymptoms]
                 .get_support_max(dt),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
                 .get_support_max(dt),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToRecovered]
                 .get_support_max(dt),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToInfectedSevere]
                 .get_support_max(dt),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToRecovered]
                 .get_support_max(dt),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSevereToInfectedCritical]
                 .get_support_max(dt),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSevereToRecovered]
                 .get_support_max(dt),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedCriticalToDead]
                 .get_support_max(dt),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedCriticalToRecovered]
                 .get_support_max(dt)});

        if (m_transitions.get_num_time_points() < (Eigen::Index)std::ceil(support_max / dt)) {
            log_error(
                "Initialization failed. Not enough time points for transitions given before start of simulation.");
        }

        parameters.check_constraints();
    }

    /**
     * @brief Calculate the number of individuals in each compartment for time 0.
     * 
     * Initial transitions are used to calculate the initial compartment sizes.
     * @param[in] dt Time discretization step size.         
     */
    void initialize(ScalarType dt);

    /**
    * @brief Computes number of Susceptibles for the current last time in m_populations.
    *
    * Number is computed using previous number of Susceptibles and the force of infection (also from previous timestep).
    * Number is stored at the matching index in m_populations.
    * @param[in] dt Time discretization step size.    
    */
    void compute_susceptibles(ScalarType dt);

    /**
     * @brief Computes size of a flow.
     * 
     * Computes size of one flow from #InfectionTransition, specified in idx_InfectionTransitions, for the current 
     * last time value in m_transitions.
     *
     * @param[in] idx_InfectionTransitions Specifies the considered flow from #InfectionTransition.
     * @param[in] idx_IncomingFlow Index of the flow in #InfectionTransition, which goes to the considered starting
     *      compartment of the flow specified in idx_InfectionTransitions. Size of considered flow is calculated via 
     *      the value of this incoming flow.
     * @param[in] dt Time step to compute flow for.
     */
    void compute_flow(int idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt);

    /**
     * @brief Sets all required flows for the current last timestep in m_transitions.
     *
     * New values are stored in m_transitions. Most values are computed via the function compute_flow().
     *
     * @param[in] dt Time step.
     */
    void flows_current_timestep(ScalarType dt);

    /**
     * @brief Computes total number of Deaths for the current last time in m_populations.
     * 
     * Number is stored in m_populations.
     *
     */
    void compute_deaths();

    /**
     * @brief Computes force of infection for the current last time in m_transitions.
     * 
     * Computed value is stored in m_forceofinfection.
     * 
     * @param[in] dt Time discretization step size.          
     * @param[in] initialization if true we are in the case of the initilization of the model. 
     *      For this we need forceofinfection at timepoint -dt which differs to usually used timepoints.
     */
    void update_forceofinfection(ScalarType dt, bool initialization = false);

    /**
     * @brief Get the size of the compartment specified in idx_InfectionState at the current last time in m_populations.
     * 
     * Calculation is reasonable for all compartments except S, R, D. 
     * Therefore, we have alternative functions for those compartments.
     *
     * @param[in] idx_InfectionState Specifies the considered #InfectionState
     * @param[in] idx_IncomingFlow Specifies the index of the infoming flow to #InfectionState in m_transitions. 
     * @param[in] idx_TransitionDistribution1 Specifies the index of the first relevant TransitionDistribution, 
     *              related to a flow from the considered #InfectionState to any other #InfectionState.
     *              This index is also used for related probability.
     * @param[in] idx_TransitionDistribution2 Specifies the index of the second relevant TransitionDistribution, 
     *              related to a flow from the considered #InfectionState to any other #InfectionState (in most cases to Recovered). 
     *              Necessary related probability is calculated via 1-probability[idx_TransitionDistribution1].
     * @param[in] dt Time discretization step size.
     */
    void compute_compartment(Eigen::Index idx_InfectionState, Eigen::Index idx_IncomingFlow,
                             int idx_TransitionDistribution1, int idx_TransitionDistribution2, ScalarType dt);

    /**
     * @brief Sets all values of remaining compartments (compartments apart from S, R, D) for the current last timestep in m_populations.
     *
     * New values are stored in m_populations. Most values are computed via the function get_size_of_compartments().
     * 
     * @param[in] dt Time discretization step size.
     */
    void other_compartments_current_timestep(ScalarType dt);

    /**
     * @brief Computes total number of Recovered for the current last time in m_populations.
     * 
     * Number is stored in m_populations.
     *
     */
    void compute_recovered();

    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    /* Attention: m_populations and m_transitions do not necessarily have the same number of time points due to the initialization part. */
    TimeSeries<ScalarType>
        m_transitions; ///< TimeSeries containing points of time and the corresponding number of transitions.
    TimeSeries<ScalarType>
        m_populations; ///< TimeSeries containing points of time and the corresponding number of people in defined #InfectionState%s.

private:
    ScalarType m_forceofinfection{0}; ///< Force of infection term needed for numerical scheme.
    ScalarType m_N{0}; ///< Total population size of the considered region.
    ScalarType m_deaths_before{0}; ///< Deaths before start of simulation (at time -m_dt).
};

} // namespace isecir
} // namespace mio

#endif // IDESECIR_MODEL_H
