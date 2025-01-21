/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/custom_index_array.h"
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
    * @param[in, out] transitions_init TimeSeries with the initial values of the number of individuals, 
    *   which transit within one timestep dt from one compartment to another.
    *   Possible transitions are specified in #InfectionTransition%s.
    *   Considered time points should have the distance dt. The last time point determines the start time t0 of the 
    *   simulation. 
    *   The time history must reach a certain point in the past so that the simulation can be performed.
    *   A warning is displayed if the condition is violated.
    * @param[in] N_init A vector, containing the populations of the considered region, for every AgeGroup.
    * @param[in] deaths_init A vector, containing the total number of deaths at time t0, for every AgeGroup.
    * @param[in] num_agegroups The number of AgeGroups.
    * @param[in] total_confirmed_cases_init A vector, containing the total confirmed cases at time t0 can be set if it 
    *   should be used for initialization, for every AgeGroup.
    */
    Model(TimeSeries<ScalarType>&& transitions_init, CustomIndexArray<ScalarType, AgeGroup> N_init,
          CustomIndexArray<ScalarType, AgeGroup> deaths_init, size_t num_agegroups,
          CustomIndexArray<ScalarType, AgeGroup> total_confirmed_cases_init = CustomIndexArray<ScalarType, AgeGroup>());

    // ---- Additional functionality such as constraint checking, setters and getters, etc. ----
    /**
    * @brief Checks constraints on model parameters and initial data.
    * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false.
    */
    bool check_constraints(ScalarType dt) const;

    /**
    * @brief Returns a flat index for the TimeSeries transitions which contains values for the #InfectionTransition%s.
    *
    * In the TimeSeries we store a vector for each time point. In this vector we store values for the different 
    * #InfectionTransition%s for every AgeGroup.
    * This function is used to get the right index in this vector for a specific AgeGroup and #InfectionTransition.
    *
    * @param[in] transition_idx Index determining which #InfectionTransition we want to evaluate.
    * @param[in] agegroup The agegroup for which we want to evaluate.
    */

    int get_transition_flat_index(Eigen::Index transition_idx, AgeGroup agegroup) const
    {
        return (static_cast<int>(size_t(agegroup)) * int(InfectionTransition::Count) + int(transition_idx));
    }

    /**
    * @brief Returns a flat index for the TimeSeries populations which contains values for the #InfectionState%s.
    *
    * In the TimeSeries we store a vector for each time point. In this vector we store values for the 
    * different #InfectionState%s for every AgeGroup.
    * This function is used to get the right index in this vector for a specific AgeGroup and #InfectionState.
    *
    * @param[in] state_idx Index at which #InfectionState we want to evaluate.
    * @param[in] agegroup The agegroup for which we want to evaluate.
    */
    int get_state_flat_index(Eigen::Index state_idx, AgeGroup agegroup) const
    {
        return (static_cast<int>(size_t(agegroup)) * int(InfectionState::Count) + int(state_idx));
    }

    /**
     * @brief Getter for the global support_max, i.e. the maximum of support_max over all TransitionDistributions.
     *
     * This determines how many initial values we need for the transitions.
     * It may be possible to run the simulation with fewer time points than the value of the global support_max, 
     * but this number ensures that it is possible.
     *
     * @param[in] dt Time step size.
     * 
     * @return Global support_max.
     */
    ScalarType get_global_support_max(ScalarType dt) const;

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
    int get_initialization_method_compartments() const
    {
        return m_initialization_method;
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

    // ---- Public parameters. ----
    ParameterSet parameters{AgeGroup(m_num_agegroups)}; ///< ParameterSet of Model Parameters.
    // Attention: populations and transitions do not necessarily have the same number of time points due to the
    // initialization part.
    TimeSeries<ScalarType>
        transitions; ///< TimeSeries containing points of time and the corresponding number of individuals transitioning from
    // one #InfectionState to another as defined in #InfectionTransition%s for every AgeGroup.
    TimeSeries<ScalarType> populations; ///< TimeSeries containing points of time and the corresponding number of
        // people in defined #InfectionState%s for every AgeGroup.
    CustomIndexArray<ScalarType, AgeGroup>
        total_confirmed_cases; ///< CustomIndexArray that contains the total number of confirmed cases at time t0 for every AgeGroup.

private:
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
     * @param[in] group The AgeGroup for which we want to compute.
     * @param[in] idx_IncomingFlow Specifies the index of the incoming flow to #InfectionState in transitions. 
     * @param[in] idx_TransitionDistribution1 Specifies the index of the first relevant TransitionDistribution, 
     *              related to a flow from the considered #InfectionState to any other #InfectionState.
     *              This index is also used for related probability.
     * @param[in] idx_TransitionDistribution2 Specifies the index of the second relevant TransitionDistribution, 
     *              related to a flow from the considered #InfectionState to any other #InfectionState (in most cases
     *               to Recovered). 
     *              Related probability is calculated via 1-probability[idx_TransitionDistribution1].
     *              Sometimes the second index is not needed, e.g., if probability[idx_TransitionDistribution1]=1.
     */

    void compute_compartment_from_flows(ScalarType dt, Eigen::Index idx_InfectionState, AgeGroup group,
                                        Eigen::Index idx_IncomingFlow, int idx_TransitionDistribution1,
                                        int idx_TransitionDistribution2 = 0);

    /**
     * @brief Computes the values of the infection compartments subset at initialization.
     *
     * The values for the compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and 
     * InfectedCritical for time t_0 are calculated using the initial data in form of flows.
     * Calculated values are stored in populations.
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
    * @brief Computes number of Susceptible%s for the current last time in populations.
    *
    * Number is computed using previous number of Susceptible%s and the force of infection (also from previous timestep).
    * Number is stored at the matching index in populations.
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
     * @param[in] group The Age group for which we want to compute the flow.
     */
    void compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt,
                      Eigen::Index current_time_index, AgeGroup group);

    /**
     * @brief Computes size of a flow for the current last time value in transitions.
     * 
     * Computes size of one flow from #InfectionTransition, specified in idx_InfectionTransitions, for the current 
     * last time value in transitions. 
     *
     * @param[in] idx_InfectionTransitions Specifies the considered flow from #InfectionTransition.
     * @param[in] idx_IncomingFlow Index of the flow in #InfectionTransition, which goes to the considered starting
     *      compartment of the flow specified in idx_InfectionTransitions. Size of considered flow is calculated via 
     *      the value of this incoming flow.
     * @param[in] dt Time step to compute flow for.
     * @param[in] group The Age Group for which we want to compute the flow.
     */
    void compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt,
                      AgeGroup group);

    /**
     * @brief Sets all required flows for the current last timestep in transitions.
     *
     * New values are stored in transitions. Most values are computed via the function compute_flow().
     *
     * @param[in] dt Time step.
     */
    void flows_current_timestep(ScalarType dt);

    /**
     * @brief Updates the values of one compartment using flows.
     *
     * New value is stored in populations. The value is calculated using the compartment size in the previous 
     * time step and the related flows of the current time step. 
     * Therefore the flows of the current time step should be calculated before using this function.
     */
    void update_compartment_from_flow(InfectionState infectionState,
                                      std::vector<InfectionTransition> const& IncomingFlows,
                                      std::vector<InfectionTransition> const& OutgoingFlows, AgeGroup group);

    /**
     * @brief Updates the values of all compartments except Susceptible at initialization.
     *
     * New values are stored in populations. The values are calculated using the compartment size in the previous 
     * time step and the related flows of the current time step. 
     * Therefore the flows of the current time step should be calculated before using this function.
     * 
     */
    void update_compartments();

    /**
     * @brief Computes force of infection for the current last time in transitions.
     * 
     * Computed value is stored in m_forceofinfection.
     * 
     * @param[in] dt Time discretization step size.          
     * @param[in] initialization If true we are in the case of the initialization of the model. 
     *      For this we need forceofinfection at time point t0-dt and not at the current last time 
     *      (given by transitions) as in the other time steps.
     */
    void compute_forceofinfection(ScalarType dt, bool initialization = false);

    // ---- Functionality to set vectors with necessary information regarding TransitionDistributions. ----
    /**
     * @brief Setter for the vector m_transitiondistributions_support_max that contains the support_max for all 
     * TransitionDistributions.
     *
     * This determines how many summands are required when calculating flows, the force of infection or compartments.
     *
     * @param[in] dt Time step size.
     */
    void set_transitiondistributions_support_max(ScalarType dt);

    /**
     * @brief Setter for the vector m_transitiondistributions_derivative that contains the approximated derivative for 
     * all TransitionDistributions for all necessary time points.
     *
     * The derivative is approximated using a backwards difference scheme.
     * The number of necessary time points for each TransitionDistribution is determined using 
     * m_transitiondistributions_support_max.
     *
     * @param[in] dt Time step size.
     */
    void set_transitiondistributions_derivative(ScalarType dt);

    /**
     * @brief Setter for the vector m_transitiondistributions_in_forceofinfection.
     *
     * When computing the force of infection, we evaluate the survival functions of the TransitionDistributions 
     * InfectedNoSymptomsToInfectedSymptoms, InfectedNoSymptomsToRecovered, InfectedSymptomsToInfectedSevere and 
     * InfectedSymptomsToRecovered, weighted by the corresponding TransitionProbabilities, at the same time points. 
     * Here, we compute these contributions to the force of infection term and store them in the vector 
     * m_transitiondistributions_in_forceofinfection so that we can access this vector for all following computations. 
     *
     * @param[in] dt Time step size.
     */
    void set_transitiondistributions_in_forceofinfection(ScalarType dt);

    // ---- Private parameters. ----
    CustomIndexArray<ScalarType, AgeGroup> m_forceofinfection; ///< Force of infection term needed for numerical scheme.
    CustomIndexArray<ScalarType, AgeGroup>
        m_N; ///< Vector containing the total population size of the considered region for every AgeGroup.
    ScalarType m_tol{1e-10}; ///< Tolerance used to calculate the maximum support of the TransitionDistributions.
    size_t m_num_agegroups; ///< Number of Age Groups.
    int m_initialization_method{0}; ///< Gives the index of the method used for the initialization of the model.
    // See also get_initialization_method_compartments() for the number code.
    CustomIndexArray<std::vector<ScalarType>, AgeGroup> m_transitiondistributions_support_max; ///<  CustomIndexArray
    // of Type AgeGroup containing a Vector containing the support_max for all TransitionDistributions.
    CustomIndexArray<std::vector<std::vector<ScalarType>>, AgeGroup>
        m_transitiondistributions_derivative; ///< CustomIndexArray
    //of Type AgeGroup containing a Vector containing the approximated derivative for all TransitionDistributions for
    // all necessary time points.
    CustomIndexArray<std::vector<std::vector<ScalarType>>, AgeGroup>
        m_transitiondistributions_in_forceofinfection; ///< CustomIndexArray
    // of Type AgeGroup containing a Vector containing the weighted TransitionDistributions for all necessary time
    // points in the force of infection term.

    // ---- Friend classes/functions. ----
    // In the Simulation class, the actual simulation is performed which is why it needs access to the here
    // defined (and private) functions to solve the model equations.
    friend class Simulation;
    // In set_initial_flows(), we compute initial flows based on RKI data using the (private) compute_flow() function
    // which is why it is defined as a friend function.
    friend IOResult<void> set_initial_flows(Model& model, ScalarType dt, std::string const& path, Date date,
                                            ScalarType scale_confirmed_cases);
};

} // namespace isecir
} // namespace mio

#endif // IDESECIR_MODEL_H
