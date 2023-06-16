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
    *   Possible transitions are specified in as InfectionTransition%s.
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
        // if (!(m_populations.get_num_time_points() > 0)) {
        //     log_error("Model construction failed. No initial time point for populations.");
        // }

        // for (int i = 0; i < (int)InfectionState::Count; i++) {
        //     if (m_populations[0][i] < 0) {
        //         log_error("Initialization failed. Initial values for populations are less than zero.");
        //     }
        // }

        if (!((int)m_transitions.get_num_elements() == (int)InfectionTransition::Count)) {
            log_error(
                "Initialization failed. Number of elements in transition vector does not match the required number.");
        }

        ScalarType global_max_support = get_global_max_support(dt);

        if (m_transitions.get_num_time_points() < (Eigen::Index)std::ceil(global_max_support / dt)) {
            log_error(
                "Initialization failed. Not enough time points for transitions given before start of simulation.");
        }

        parameters.check_constraints();
    }

    // define function that computes flows needed for initalization of IDE for a given result/compartments of the ODE model
    // we assume that the ODE simulation starts at t0=0
    void compute_initial_flows_from_compartments(mio::TimeSeries<ScalarType> secihurd_ode, ScalarType t0_ide,
                                                 ScalarType dt)
    {
        int num_transitions = (int)mio::isecir::InfectionTransition::Count;

        // get (global) max_support to determine how many flows in the past we have to compute
        ScalarType global_max_support         = this->get_global_max_support(dt);
        Eigen::Index global_max_support_index = std::ceil(global_max_support / dt);
        std::cout << "Global max_support: " << global_max_support << "\n";

        // remove time point
        this->m_transitions.remove_last_time_point();

        ScalarType t0_ide_index = std::ceil(t0_ide / dt);
        unused(secihurd_ode);
        // flow from S to E for -6*global_max_support, ..., 0 (directly from compartments)
        // add time points to init_transitions here
        for (int i = t0_ide_index - 6 * global_max_support_index + 1; i <= t0_ide_index; i++) {
            this->m_transitions.add_time_point(i, mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0));
            // std::cout << "i: " << i << "\n";
            this->m_transitions.get_last_value()[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)] =
                secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::Susceptible)] -
                secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::Susceptible)];
        }

        // then use compute_flow function to compute following flows

        Eigen::Index start_shift = t0_ide_index - 6 * global_max_support_index;

        // flow from E to C for -5*global_max_support, ..., 0
        for (int i = t0_ide_index - 5 * global_max_support_index + 1; i <= t0_ide_index; i++) {
            this->compute_flow(1, 0, dt, true, i - start_shift);
        }

        // flow from C to I for -4*global_max_support, ..., 0
        for (int i = t0_ide_index - 4 * global_max_support_index + 1; i <= t0_ide_index; i++) {
            // C to I
            this->compute_flow(2, 1, dt, true, i - start_shift);
        }

        // flow from I to H for -3*global_max_support, ..., 0
        for (int i = t0_ide_index - 3 * global_max_support_index + 1; i <= t0_ide_index; i++) {
            // I to H
            this->compute_flow(4, 2, dt, true, i - start_shift);
        }

        // flow from H to U for -2*global_max_support, ..., 0
        for (int i = t0_ide_index - 2 * global_max_support_index + 1; i <= t0_ide_index; i++) {
            // H to U
            this->compute_flow(6, 4, dt, true, i - start_shift);
        }

        // flow from U to D and C, I, H, U to R for -1*global_max_support, ..., 0
        for (int i = t0_ide_index - 1 * global_max_support_index + 1; i <= t0_ide_index; i++) {
            // U to D
            this->compute_flow(8, 6, dt, true, i - start_shift);
            // C to R
            this->compute_flow(3, 1, dt, true, i - start_shift);
            // I to R
            this->compute_flow(5, 2, dt, true, i - start_shift);
            // H to R
            this->compute_flow(7, 4, dt, true, i - start_shift);
            // U to R
            this->compute_flow(9, 6, dt, true, i - start_shift);
        }
    }

    /**
     * @brief Calculate the number of individuals in each compartment for time 0.
     * 
     * Initial transitions are used to calculate the initial compartment sizes.
     * @param[in] dt Time discretization step size.         
     */
    void initialize_solver(ScalarType dt);

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
     * Computes size of one flow from InfectionTransition, specified in idx_InfectionTransitions, for the current 
     * last time value in m_transitions.
     *
     * @param[in] idx_InfectionTransitions Specifies the considered flow from InfectionTransition.
     * @param[in] idx_IncomingFlow Index of the flow in InfectionTransition, which goes to the considered starting
     *      compartment of the flow specified in idx_InfectionTransitions. Size of considered flow is calculated via 
     *      the value of this incoming flow.
     * @param[in] dt Time step to compute flow for.
     */
    void compute_flow(int idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt,
                      bool initial_flow = false, Eigen::Index current_time_index = 0);

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
     * @param[in] idx_InfectionState Specifies the considered InfectionState
     * @param[in] idx_IncomingFlow Specifies the index of the infoming flow to InfectionState in m_transitions. 
     * @param[in] idx_TransitionDistribution1 Specifies the index of the first relevant TransitionDistribution, 
     *              related to a flow from the considered InfectionState to any other InfectionState.
     *              This index is also used for related probability.
     * @param[in] idx_TransitionDistribution2 Specifies the index of the second relevant TransitionDistribution, 
     *              related to a flow from the considered InfectionState to any other InfectionState (in most cases to Recovered). 
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

    /**
     * @brief Getter for the global max_support, i.e. the maximum of max_support over all TransitionDistributions.
     *
     * This determines how many inital values we need for the flows.
     *
     * @param[in] dt Time step size.
     * 
     * @return Global max_support.
     *
     */
    ScalarType get_global_max_support(ScalarType dt) const;

    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    /* Attention: m_populations and m_transitions do not necessarily have the same number of time points due to the initialization part. */
    TimeSeries<ScalarType>
        m_transitions; ///< TimeSeries containing points of time and the corresponding number of transitions.
    TimeSeries<ScalarType>
        m_populations; ///< TimeSeries containing points of time and the corresponding number of people in defined InfectionState%s.

private:
    ScalarType m_forceofinfection{0}; ///< Force of infection term needed for numerical scheme.
    ScalarType m_N{0}; ///< Total population size of the considered region.
    ScalarType m_deaths_before{0}; ///< Deaths before start of simulation (at time -m_dt).
};

} // namespace isecir
} // namespace mio

#endif // IDESECIR_MODEL_H
