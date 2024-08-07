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

#include "abm/virus_variant.h"
#include "ide_secir/parameters.h"
#include "ide_secir/infection_state.h"
#include "ode_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"

#include "vector"

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
    *   Possible transitions are specified in #InfectionTransition%s.
    *   Considered time points should have the distance dt. The last time point determines the start time t0 of the 
    *   simulation. 
    *   The time history must reach a certain point in the past so that the simulation can be performed.
    *   A warning is displayed if the condition is violated.
    * @param[in] N_init The population of the considered region.
    * @param[in] deaths The total number of deaths at time t0.
    * @param[in] total_confirmed_cases Total confirmed cases at time t0 can be set if it should be used for initialization.
    * @param[in, out] Parameterset_init Used Parameters for simulation. 
    */

    Model(TimeSeries<ScalarType>&& init, ScalarType N_init, ScalarType deaths, ScalarType total_confirmed_cases = 0,
          const ParameterSet& Parameterset_init = ParameterSet());

    // ---- Functionality to calculate the sizes of the compartments for time t0. ----
    /**
     * @brief Compute the compartment specified in idx_InfectionState at the current time -- only using historic flow 
     * values and disrespecting potential, previous compartment value.
     * 
     * The computation is meaningful for all compartments except Susceptible, Recovered and #Death
     * and mostly needed for initialization. 
     * For Susceptible, Recovered and Dead, use corresponding alternative functions.
     *
     * @param[in] dt Time discretization step size.
     * @param[in] idx_InfectionState Specifies the considered #InfectionState
     * @param[in] idx_IncomingFlow Specifies the index of the infoming flow to #InfectionState in m_transitions. 
     * @param[in] idx_TransitionDistribution1 Specifies the index of the first relevant TransitionDistribution, 
     *              related to a flow from the considered #InfectionState to any other #InfectionState.
     *              This index is also used for related probability.
     * @param[in] idx_TransitionDistribution2 Specifies the index of the second relevant TransitionDistribution, 
     *              related to a flow from the considered #InfectionState to any other #InfectionState (in most cases
     *               to Recovered). 
     *              Related probability is calculated via 1-probability[idx_TransitionDistribution1].
     *              Sometimes the second index is not needed, e.g., if probability[idx_TransitionDistribution1]=1.
     */
    void compute_compartment_from_flows(ScalarType dt, Eigen::Index idx_InfectionState, Eigen::Index idx_IncomingFlow,
                                        int idx_TransitionDistribution1, int idx_TransitionDistribution2 = 0);

    /**
     * @brief Computes the values of the infection compartments subset at initialization.
     *
     * The values for the compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical
     * for time t_0 are calculated using the initial data in form of flows.
     * Calculated values are stored in m_populations.
     * 
     * @param[in] dt Time discretization step size.
     */
    void initial_compute_compartments_infection(ScalarType dt);

    /**
     * @brief Computes the values of compartments at initialization.
     * 
     * The initialization method is selected automatically based on the different values that need to be set beforehand.
     * Infection compartments are always computed through historic flows.
     * Initialization methods for Susceptible and Recovered are tested in the following order:
     * 1.) If a positive number for the total number of confirmed cases is set, Recovered is set according to that 
     * value and Susceptible%s are derived.
     * 2.) If Susceptible%s are set, Recovered will be derived.
     * 3.) If Recovered are set directly, Susceptible%s are derived.
     * 4.) If none of the above is set with positive value, the force of infection is used as in Messina et al (2021) 
     * to set the Susceptible%s.
     *
     * The function calculate_compartments_initialization() is used in every method for the compartments Exposed, 
     * InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.
     * In addition, the force of infection is calculated for start time t_0.
     *
     * @param[in] dt Time discretization step size.      
     */
    void initial_compute_compartments(ScalarType dt);

    // ---- Functionality for the iterations of a simulation. ----
    /**
    * @brief Computes number of Susceptible%s for the current last time in m_populations.
    *
    * Number is computed using previous number of Susceptible%s and the force of infection (also from previous timestep).
    * Number is stored at the matching index in m_populations.
    * @param[in] dt Time discretization step size.    
    */
    void compute_susceptibles(ScalarType dt);

    /**
     * @brief Computes size of a flow.
     * 
     * Computes size of one flow from #InfectionTransition, specified in idx_InfectionTransitions, for the time 
     * index current_time_index.
     *
     * @param[in] idx_InfectionTransitions Specifies the considered flow from #InfectionTransition.
     * @param[in] idx_IncomingFlow Index of the flow in #InfectionTransition, which goes to the considered starting
     *      compartment of the flow specified in idx_InfectionTransitions. Size of considered flow is calculated via 
     *      the value of this incoming flow.
     * @param[in] dt Time step to compute flow for.
     * @param[in] current_time_index The time index the flow should be computed for.
     */
    void compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt,
                      Eigen::Index current_time_index);

    /**
     * @brief Computes size of a flow for the current last time value in m_transitions.
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
    void compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt);

    /**
     * @brief Sets all required flows for the current last timestep in m_transitions.
     *
     * New values are stored in m_transitions. Most values are computed via the function compute_flow().
     *
     * @param[in] dt Time step.
     */
    void flows_current_timestep(ScalarType dt);

    /**
     * @brief Updates the values of all compartments except Susceptible at initialization.
     *
     * New values are stored in m_populations. The values are calculated using the compartment size in the previous 
     * time step and the related flows of the current time step. 
     * Therefore the flows of the current time step should be calculated before using this function.
     * 
     */
    void update_compartments();

    /**
     * @brief Computes force of infection for the current last time in m_transitions.
     * 
     * Computed value is stored in m_forceofinfection.
     * 
     * @param[in] dt Time discretization step size.          
     * @param[in] initialization If true we are in the case of the initialization of the model. 
     *      For this we need forceofinfection at time point t0-dt and not at the current last time 
     *      (given by m_transitions) as in the other time steps.
     */
    void compute_forceofinfection(ScalarType dt, bool initialization = false);

    // ---- Additional functionality such as constraint checking, setters and getters, etc. ----
    /**
    * @brief Checks constraints on model parameters and initial data.
    * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false.
    */
    bool check_constraints(ScalarType dt) const
    {
        if (!((int)m_transitions.get_num_elements() == (int)InfectionTransition::Count)) {
            log_error("A variable given for model construction is not valid. Number of elements in transition vector "
                      "does not match the required number.");
            return true;
        }

        for (int i = 0; i < (int)InfectionState::Count; i++) {
            if (m_populations[0][i] < 0) {
                log_error("Initialization failed. Initial values for populations are less than zero.");
                return true;
            }
        }

        // It may be possible to run the simulation with fewer time points, but this number ensures that it is possible.
        if (m_transitions.get_num_time_points() < (Eigen::Index)std::ceil(get_global_support_max(dt) / dt)) {
            log_error("Initialization failed. Not enough time points for transitions given before start of simulation. "
                      "{} time points are required.",
                      (Eigen::Index)std::ceil(get_global_support_max(dt) / dt));
            return true;
        }

        for (int i = 0; i < m_transitions.get_num_time_points(); i++) {
            for (int j = 0; j < (int)InfectionTransition::Count; j++) {
                if (m_transitions[i][j] < 0) {
                    log_error("Initialization failed. One or more initial value for transitions is less than zero.");
                    return true;
                }
            }
        }
        if (m_transitions.get_last_time() != m_populations.get_last_time()) {
            log_error(
                "Last time point of TimeSeries for transitions {} does not match last time point of TimeSeries for "
                "compartments {}. For a meaningful simulation, these two points must match.",
                m_transitions.get_last_time(), m_populations.get_last_time());
            return true;
        }

        if (m_populations.get_num_time_points() != 1) {
            log_error("The TimeSeries for the compartments contains more than one time point. It is unclear how to "
                      "initialize.");
            return true;
        }

        return parameters.check_constraints();
    }

    /**
     * @brief Getter for the global support_max, i.e. the maximum of support_max over all TransitionDistributions.
     *
     * This determines how many inital values we need for the flows.
     * It may be possible to run the simulation with fewer time points than the value of the global support_max, 
     * but this number ensures that it is possible.
     *
     * @param[in] dt Time step size.
     * 
     * @return Global support_max.
     *
     */
    ScalarType get_global_support_max(ScalarType dt) const;

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
    * @brief Returns the index of the automatically selected initialization method.
    *
    * The initialization method is selected automatically based on the different values that need to be set beforehand.
    * Infection compartments are always computed through historic flow.
    * Initialization methods for Susceptible and Recovered are tested in the following order:
    * 1.) If a positive number for the total number of confirmed cases is set, Recovered is set according to that value
    *    and Susceptible%s are derived.
    * 2.) If Susceptible%s are set, Recovered will be derived.
    * 3.) If Recovered are set directly, Susceptible%s are derived.
    * 4.) If none of the above is set with positive value, the force of infection is used as in Messina et al (2021) to
    *     set the Susceptible%s.
    *
    * @return Index representing the initialization method.
    */
    int get_initialization_method_compartments()
    {
        return m_initialization_method;
    }

    // ---- Public parameters. ----
    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    // Attention: m_populations and m_transitions do not necessarily have the same number of time points due to the
    // initialization part.
    TimeSeries<ScalarType>
        m_transitions; ///< TimeSeries containing points of time and the corresponding number of transitions.
    TimeSeries<ScalarType> m_populations; ///< TimeSeries containing points of time and the corresponding number of
        // people in defined #InfectionState%s.
    ScalarType m_total_confirmed_cases{0}; ///< Total number of confirmed cases at time t0.

private:
    /**
     * @brief Updates the values of one compartment using flows.
     *
     * New value is stored in m_populations. The value is calculated using the compartment size in the previous 
     * time step and the related flows of the current time step. 
     * Therefore the flows of the current time step should be calculated before using this function.
     */
    void update_compartment_from_flow(InfectionState infectionState,
                                      std::vector<InfectionTransition> const& IncomingFlows,
                                      std::vector<InfectionTransition> const& OutgoingFlows);

    // ---- Private parameters. ----
    ScalarType m_forceofinfection{0}; ///< Force of infection term needed for numerical scheme.
    ScalarType m_N{0}; ///< Total population size of the considered region.
    ScalarType m_tol{1e-10}; ///< Tolerance used to calculate the maximum support of the TransitionDistributions.
    int m_initialization_method{0}; ///< Gives the index of the method used for the initialization of the model.
        // See also get_initialization_method_compartments() for the number code.
};

} // namespace isecir
} // namespace mio

#endif // IDESECIR_MODEL_H
