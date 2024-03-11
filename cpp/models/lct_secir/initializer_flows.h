/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Lena Ploetzke
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

#ifndef LCTSECIR_INITIALIZER_H
#define LCTSECIR_INITIALIZER_H

#include "memilio/config.h"
#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/logging.h"
#include "memilio/epidemiology/state_age_function.h"

namespace mio
{
namespace lsecir
{

/**
 * @brief Class that can be used to compute an initialization vector out of flows for a LCT Model.
 *
 * @tparam Model is expected to be an LCT-SECIR model defined in models/lct_secir/model.h.
 */
template <typename Model>
class Initializer
{
public:
    using LctState = typename Model::LctState;

    /**
     * @brief Constructs a new Initializer object.
     *
     * @param[in] flows Initalizing TimeSeries with flows fitting to these defined in InfectionTransition. 
     *      Timesteps should be equidistant.
     * @param[in, out] model The LCT-SECIR model for which the initialization should be performed.
     */
    Initializer(TimeSeries<ScalarType>&& flows, Model& model)
        : m_flows(std::move(flows))
        , m_model(model)
    {
        m_dt = m_flows.get_time(1) - m_flows.get_time(0);
    }

    /**
     * @brief Core function of Initializer.
     *
     * Computes a vector that can be used for the initalization of an LCT model with the number of persons for each subcompartment.
     * The initial value vector is updated in the model.
     *
     * @param[in] total_population The total size of the considered population.
     * @param[in] deaths Number of deceased people from the disease at time 0.
     * @param[in] total_confirmed_cases Total number of confirmed cases at time 0.
     * @return Returns true if one (or more) constraint(s) of the model or the initial flows are not satisfied, otherwise false. 
     */
    bool compute_initialization_vector(ScalarType total_population, ScalarType deaths, ScalarType total_confirmed_cases)
    {

        Eigen::VectorXd init(LctState::Count);

        //E
        init.segment(LctState::template get_first_index<InfectionState::Exposed>(),
                     LctState::template get_num_subcompartments<InfectionState::Exposed>()) =
            compute_compartment<InfectionState::Exposed>(Eigen::Index(InfectionTransition::SusceptibleToExposed),
                                                         1 / m_model.parameters.template get<TimeExposed>());
        //C
        init.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>()) =
            compute_compartment<InfectionState::InfectedNoSymptoms>(
                Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                1 / m_model.parameters.template get<TimeInfectedNoSymptoms>());
        //I
        init.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>()) =
            compute_compartment<InfectionState::InfectedSymptoms>(
                Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                1 / m_model.parameters.template get<TimeInfectedSymptoms>());
        //H
        init.segment(LctState::template get_first_index<InfectionState::InfectedSevere>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedSevere>()) =
            compute_compartment<InfectionState::InfectedSevere>(
                Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                1 / m_model.parameters.template get<TimeInfectedSevere>());
        //U
        init.segment(LctState::template get_first_index<InfectionState::InfectedCritical>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedCritical>()) =
            compute_compartment<InfectionState::InfectedCritical>(
                Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                1 / m_model.parameters.template get<TimeInfectedCritical>());
        //R
        // Number of recovered is equal to the cumulative number of confirmed cases minus the number of people who are infected at the moment.
        init[LctState::Count - 2] =
            total_confirmed_cases -
            init.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                         LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>() +
                             LctState::template get_num_subcompartments<InfectionState::InfectedSevere>() +
                             LctState::template get_num_subcompartments<InfectionState::InfectedCritical>())
                .sum() -
            deaths;

        //S
        init[0] = total_population - init.segment(1, LctState::Count - 2).sum() - deaths;

        //D
        init[LctState::Count - 1] = deaths;

        for (int i = 0; i < (int)LctState::Count; i++) {
            if (init[i] < 0) {
                log_error("Initialization failed. Not enough time points for the transitions are given. More than {} "
                          "are needed but "
                          "just {} are given.",
                          -init[i], m_flows.get_num_time_points());
                return true;
            }
        }

        // Update initial value vector of the mpdel.
        m_model.set_initial_values(init);

        // Check if constraints were fulfilled.
        return check_constraints();
    }

    /**
     * @brief Setter for the tolerance used to calculate the maximum support of ErlangDensity%s.
     *
     * @param[in] new_tol New tolerance.
     */
    void set_tol_for_support_max(ScalarType new_tol)
    {
        m_tol = new_tol;
    }

private:
    /**
     * @brief Checks constraints of the Initializer inclusive check for the model.
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false. 
     */
    bool check_constraints() const
    {
        if (!((int)InfectionTransition::Count == (int)m_flows.get_num_elements())) {
            log_error("Initial condition size does not match subcompartments.");
            return true;
        }

        if (!(floating_point_equal(0., m_flows.get_last_time(), 1e-8, 1e-14))) {
            log_error("Last time point in flows has to be 0.");
            return true;
        }

        for (int i = 1; i < m_flows.get_num_time_points(); i++) {
            if (!(floating_point_equal(m_dt, m_flows.get_time(i) - m_flows.get_time(i - 1), 1e-8, 1e-14))) {
                log_error("Time points in flows have to be equidistant.");
                return true;
            }
        }

        if (!(m_dt < 1)) {
            log_warning("Step size was set very large. The result could be distorted.");
            return true;
        }

        // Check if model is valid and the calculated initial value vector was valid.
        return m_model.check_constraints();
    }

    /**
     * @brief Computes a vector with the number of people in each subcompartment for one InfectionState.
     *
     * With this function, partial result of compute_initialization_vector are achieved.
     *
     * @param[in] idx_incoming_flow Index of the flow which is relevant for the calculation, so the flow to the InfectionState.
     * @param[in] transition_rate Specifies the transition rate of the InfectionState. Is equal to 1 / (expected Time in the InfectionState).
     * @tparam State The InfectionState for which the partial result should be calculated.
     * @return Vector with a possible initialization for the subcompartments of the InfectionState State. 
     *      Subcompartment is set to -1 if calculation was not possible.
     */
    template <InfectionState State>
    Eigen::VectorXd compute_compartment(Eigen::Index idx_incoming_flow, ScalarType transition_rate) const
    {
        int num_infectionstates = LctState::template get_num_subcompartments<State>();
        Eigen::VectorXd subcompartments(num_infectionstates);
        // Initialize relevant density for the InfectionState.
        // For the first subcompartment a shape parameter of one is needed.
        ErlangDensity erlang(1, 1. / (num_infectionstates * transition_rate));

        // Initialize other relevant parameters.
        ScalarType calc_time{0};
        Eigen::Index calc_time_index{0};
        ScalarType state_age{0};
        ScalarType sum{0};
        Eigen::Index num_time_points = m_flows.get_num_time_points();

        // Calculate number of people in each subcomaprtment.
        for (int j = 0; j < num_infectionstates; j++) {
            // For subcompartment number j+1, shape parameter j+1 is needed.
            erlang.set_parameter(j + 1);

            // Determine relevant calculation area and corresponding index.
            calc_time       = erlang.get_support_max(m_dt, m_tol);
            calc_time_index = (Eigen::Index)std::ceil(calc_time / m_dt) - 1;

            if (num_time_points < calc_time_index) {
                // Initialization failed. Not enough time points for the transitions are given.
                // Store the number of necessary time points in the vector and throw one error at the compute_initialization_vector() function,
                // so that the error is not displayed several times.
                subcompartments[j] = -calc_time_index;
            }
            else {

                // Approximate integral with non-standard scheme.
                for (Eigen::Index i = num_time_points - calc_time_index; i < num_time_points; i++) {
                    state_age = (num_time_points - i) * m_dt;
                    sum += erlang.eval(state_age) * m_flows[i][idx_incoming_flow];
                }
                subcompartments[j] = 1 / (num_infectionstates * transition_rate) * sum;
            }
            sum = 0;
        }
        return subcompartments;
    }

    TimeSeries<ScalarType> m_flows; ///< TimeSeries with the flows which are used to calculate the initial vector.
    Model& m_model; ///< The LCT-SECIR model for which the initialization should be performed.
    ScalarType m_dt{}; ///< Step size of the times in m_flows and time step for the approximation of the integral.
    ScalarType m_tol{1e-10}; ///< Tolerance used to calculate the maximum support of the ErlangDensity%s.
};

} // namespace lsecir
} // namespace mio

#endif