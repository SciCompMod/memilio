/* 
* Copyright (C) 2020-2024 MEmilio
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
    * @brief Constructor to create an IDE-SECIR model.
    *
    * @param[in, out] init TimeSeries with the initial values of the number of individuals, 
    *   which transit within one timestep dt from one compartment to another.
    *   Possible transitions are specified in as #InfectionTransition%s.
    *   Considered points of times should have the distance dt and the last time point should be 0. 
    *   The time history must reach a certain point in the past so that the simulation can be performed.
    *   A warning is displayed if the condition is violated.
    * @param[in] N_init The population of the considered region.
    * @param[in] deaths The total number of deaths at the time zero.
    * @param[in] total_confirmed_cases Total confirmed cases at time t0 can be set if it should be used for initialisation.
    * @param[in, out] Parameterset_init Used Parameters for simulation. 
    */
    Model(TimeSeries<ScalarType>&& init, ScalarType N_init, ScalarType deaths, ScalarType total_confirmed_cases = 0,
          const ParameterSet& Parameterset_init = ParameterSet());

    // ---- Functionality to calculate the sizes of the compartments for time zero. ----
    /**
     * @brief Get the size of the compartment specified in idx_InfectionState at the current last time in m_populations.
     * 
     * Calculation is reasonable for all compartments except S, R, D.
     * Function is currently just used to calculate the compartment sizes at time t_0=0. It could also be used for other time steps.
     *
     * @param[in] idx_InfectionState Specifies the considered #InfectionState
     * @param[in] idx_IncomingFlow Specifies the index of the infoming flow to #InfectionState in m_transitions. 
     * @param[in] idx_TransitionDistribution1 Specifies the index of the first relevant TransitionDistribution, 
     *              related to a flow from the considered #InfectionState to any other #InfectionState.
     *              This index is also used for related probability.
     * @param[in] idx_TransitionDistribution2 Specifies the index of the second relevant TransitionDistribution, 
     *              related to a flow from the considered #InfectionState to any other #InfectionState (in most cases to Recovered). 
     *              Necessary related probability is calculated via 1-probability[idx_TransitionDistribution1].
     *              If the second index is not needed, eg if probability[idx_TransitionDistribution1]=1, 
     *              just use an arbitrary legal index.
     * @param[in] dt Time discretization step size.
     */
    void compute_compartment_from_flows(Eigen::Index idx_InfectionState, Eigen::Index idx_IncomingFlow,
                                        int idx_TransitionDistribution1, int idx_TransitionDistribution2,
                                        ScalarType dt);

    /**
     * @brief Sets the values of the compartments listed below for time t_0=0 in m_populations during initialization.
     *
     * The values for the compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical
     * for time t_0=0 are calculated using the initial data in form of flows.
     * Calculated values are stored in m_populations. The values are computed via the function compute_compartment_from_flows().
     * 
     * @param[in] dt Time discretization step size.
     */
    void compute_compartments_initialization(ScalarType dt);

    /**
     * @brief The initial values for the flows and additional information is used to calculate the compartment sizes at time t_0=0.
     * 
     * There are different methods to calculate the compartment sizes of the compartments Susceptibles and Recovered depending on the information given. 
     * See also the function get_initialization_method_compartments().
     * The function calculate_compartments_initialization() is used in every method for the compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.
     * In addition, the force of infection is calculated for time t_0=0.
     * 
     * @param[in] dt Time discretization step size.         
     */
    void calculate_initial_compartment_sizes(ScalarType dt);

    // ---- Functionality for the iterations of a simulation. ----
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
     * @brief Sets the values of compartments (all apart from S) for the current last timestep in m_populations.
     *
     * New values are stored in m_populations. The values are calculated using the compartment size in the previous time step and the related flows of the current time step. 
     * Therefore the flows of the current time step should be calculated before using this function.
     * 
     */
    void update_compartments_current_timestep();

    /**
     * @brief Computes force of infection for the current last time in m_transitions.
     * 
     * Computed value is stored in m_forceofinfection.
     * 
     * @param[in] dt Time discretization step size.          
     * @param[in] initialization if true we are in the case of the initilization of the model. 
     *      For this we need forceofinfection at timepoint -dt which differs to usually used timepoints.
     */
    void compute_forceofinfection(ScalarType dt, bool initialization = false);

    // ---- Additional functionality such as constraint checking, setters and getters, etc. ----
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
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToRecovered]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToInfectedSevere]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToRecovered]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSevereToInfectedCritical]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSevereToRecovered]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedCriticalToDead]
                 .get_support_max(dt, m_tol),
             parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedCriticalToRecovered]
                 .get_support_max(dt, m_tol)});

        if (m_transitions.get_num_time_points() < (Eigen::Index)std::ceil(support_max / dt)) {
            log_error(
                "Initialization failed. Not enough time points for transitions given before start of simulation.");
        }

        parameters.check_constraints();
    }

    /**
     * @brief Setter for the tolerance used to calculate the maximum support of the TransitionDistributions.
     *
     * @param[in] new_tol New tolerance.
     */
    void set_tol_for_support_max(ScalarType new_tol)
    {
        m_tol = new_tol;
    }

    /**
     * @brief Specifies a number associated with the method used for initialization.
     *
     * @returns 0 if the initialization method has not yet been selected,
     *      1 if the method using the total number of confirmed cases at time 0 is used,
     *      2 if the force of infection method is used,
     *      3 if the initialization is calculated using a prior set value for S,
     *      4 if the initialization is calculated using a prior set value for R and
     *      -1 if the initialization was not possible using any of the methods.
     */
    int get_initialization_method_compartments()
    {
        return m_initialization_method;
    }

    // ---- Public parameters. ----
    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    /* Attention: m_populations and m_transitions do not necessarily have the same number of time points due to the initialization part. */
    TimeSeries<ScalarType>
        m_transitions; ///< TimeSeries containing points of time and the corresponding number of transitions.
    TimeSeries<ScalarType>
        m_populations; ///< TimeSeries containing points of time and the corresponding number of people in defined #InfectionState%s.

private:
    // ---- Private parameters. ----
    ScalarType m_forceofinfection{0}; ///< Force of infection term needed for numerical scheme.
    ScalarType m_N{0}; ///< Total population size of the considered region.
    ScalarType m_total_confirmed_cases{0}; ///< Total number of confirmed cases at time t0.
    ScalarType m_tol{1e-10}; ///< Tolerance used to calculate the maximum support of the TransitionDistributions.
    int m_initialization_method{
        0}; ///< Gives the index of the method used for the initialization of the model. See also get_initialization_method_compartments() for the number code.
};

} // namespace isecir
} // namespace mio

#endif // IDESECIR_MODEL_H
