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
 *  The initialization method is based on the fact that the LCT model is a special case of an IDE model with special 
 *  choices of stay time distributions (so-called Erlang distributions). Accordingly, the method for calculating 
 *  initial values for the compartments from given flows/transitions is taken from the IDE-SECIR model.
 *  We cannot use the functionality of the IDE model directly, as we have to calculate a division into sub-compartments
 *  and not only the sizes of the compartments.
*   See also the IDE-SECIR model for the general method and for a better understanding of flows/transitions.
 *
 * @tparam Model is expected to be an LCT-SECIR model defined in models/lct_secir/model.h.
 */
template <typename Model>
class Initializer
{
public:
    using LctState       = typename Model::LctState;
    using InfectionState = typename LctState::InfectionState;

    /**
     * @brief Constructs a new Initializer object.
     *
     * @param[in] flows Initalizing TimeSeries with flows fitting to these defined in InfectionTransition. 
     *      Timesteps should be equidistant and the values should be non-negative. 
     *      The time history has to be long enough so that it is possible to calculate the initial vector. 
     *      The length of the required time history depends on the Erlang densities used to compute the initial vector.
     *      An error is thrown if this condition is violated.
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

        // Exposed.
        init.segment(LctState::template get_first_index<InfectionState::Exposed>(),
                     LctState::template get_num_subcompartments<InfectionState::Exposed>()) =
            compute_compartment<InfectionState::Exposed>((Eigen::Index)InfectionTransition::SusceptibleToExposed,
                                                         1 / m_model.parameters.template get<TimeExposed>());

        if ((init.segment(LctState::template get_first_index<InfectionState::Exposed>(),
                          LctState::template get_num_subcompartments<InfectionState::Exposed>())
                 .array() < 0)
                .any()) {
            // Something went wrong while calculating the initial values for the subcompartments.
            // See the compute_compartment() function for details.
            return true;
        }
        // InfectedNoSymptoms.
        init.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>()) =
            compute_compartment<InfectionState::InfectedNoSymptoms>(
                (Eigen::Index)InfectionTransition::ExposedToInfectedNoSymptoms,
                1 / m_model.parameters.template get<TimeInfectedNoSymptoms>());

        if ((init.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                          LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
                 .array() < 0)
                .any()) {
            // Something went wrong while calculating the initial values for the subcompartments.
            // See the compute_compartment() function for details.
            return true;
        }
        // InfectedSymptoms.
        init.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>()) =
            compute_compartment<InfectionState::InfectedSymptoms>(
                (Eigen::Index)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
                1 / m_model.parameters.template get<TimeInfectedSymptoms>());

        if ((init.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                          LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
                 .array() < 0)
                .any()) {
            // Something went wrong while calculating the initial values for the subcompartments.
            // See the compute_compartment() function for details.
            return true;
        }
        // InfectedSevere.
        init.segment(LctState::template get_first_index<InfectionState::InfectedSevere>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedSevere>()) =
            compute_compartment<InfectionState::InfectedSevere>(
                (Eigen::Index)InfectionTransition::InfectedSymptomsToInfectedSevere,
                1 / m_model.parameters.template get<TimeInfectedSevere>());

        if ((init.segment(LctState::template get_first_index<InfectionState::InfectedSevere>(),
                          LctState::template get_num_subcompartments<InfectionState::InfectedSevere>())
                 .array() < 0)
                .any()) {
            // Something went wrong while calculating the initial values for the subcompartments.
            // See the compute_compartment() function for details.
            return true;
        }
        // InfectedCritical.
        init.segment(LctState::template get_first_index<InfectionState::InfectedCritical>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedCritical>()) =
            compute_compartment<InfectionState::InfectedCritical>(
                (Eigen::Index)InfectionTransition::InfectedSevereToInfectedCritical,
                1 / m_model.parameters.template get<TimeInfectedCritical>());

        if ((init.segment(LctState::template get_first_index<InfectionState::InfectedCritical>(),
                          LctState::template get_num_subcompartments<InfectionState::InfectedCritical>())
                 .array() < 0)
                .any()) {
            // Something went wrong while calculating the initial values for the subcompartments.
            // See the compute_compartment() function for details.
            return true;
        }
        // Recovered.
        // Number of recovered is equal to the cumulative number of confirmed cases minus the number of people who are infected at the moment.
        init[LctState::Count - 2] =
            total_confirmed_cases -
            init.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                         LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>() +
                             LctState::template get_num_subcompartments<InfectionState::InfectedSevere>() +
                             LctState::template get_num_subcompartments<InfectionState::InfectedCritical>())
                .sum() -
            deaths;

        // Susceptibles.
        init[0] = total_population - init.segment(1, LctState::Count - 2).sum() - deaths;

        // Dead.
        init[LctState::Count - 1] = deaths;

        for (size_t i : {(size_t)0, LctState::Count - 2, LctState::Count - 1}) {
            if (init[i] < 0) {
                log_error("Initialization failed. Values of total_confirmed_cases, deaths and total_population do not "
                          "match the specified values for the transitions, leading to at least one negative result.");
                return true;
            }
        }

        // Update initial value vector of the model.
        for (size_t i = 0; i < LctState::Count; i++) {
            m_model.populations[mio::Index<LctState>(i)] = init[i];
        }

        // Check if constraints are fulfilled.
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
     * @brief Checks constraints of the Initializer including checks for the model.
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false. 
     */
    bool check_constraints() const
    {
        if (!((Eigen::Index)InfectionTransition::Count == m_flows.get_num_elements())) {
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

        // Check if model is valid and the calculated initial value vector is valid.
        m_model.check_constraints();
        return false;
    }

    /**
     * @brief Computes a vector with the number of people in each subcompartment for one InfectionState.
     *
     * This function is used in compute_initialization_vector.
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
        size_t num_subcompartments = LctState::template get_num_subcompartments<State>();
        Eigen::VectorXd subcompartments(num_subcompartments);
        // Initialize relevant density for the InfectionState.
        // For the first subcompartment a shape parameter of one is needed.
        ErlangDensity erlang(1, 1. / (num_subcompartments * transition_rate));

        // Initialize other relevant parameters.
        ScalarType calc_time{0};
        Eigen::Index calc_time_index{0};
        ScalarType state_age{0};
        ScalarType sum{0};
        Eigen::Index num_time_points = m_flows.get_num_time_points();

        // Calculate number of people in each subcompartment.
        for (size_t j = 0; j < num_subcompartments; j++) {
            // For subcompartment number j+1, shape parameter j+1 is needed.
            erlang.set_distribution_parameter((ScalarType)j + 1);

            // Determine relevant calculation area and corresponding index.
            calc_time       = erlang.get_support_max(m_dt, m_tol);
            calc_time_index = (Eigen::Index)std::ceil(calc_time / m_dt) - 1;

            if (num_time_points < calc_time_index) {
                log_error("Initialization failed. Not enough time points for the transitions are given. More than {} "
                          "are needed but just {} are given.",
                          calc_time_index, num_time_points);
                subcompartments[j] = -1.;
                return subcompartments;
            }
            else {
                // Approximate integral with non-standard scheme.
                for (Eigen::Index i = num_time_points - calc_time_index; i < num_time_points; i++) {
                    state_age = (num_time_points - i) * m_dt;
                    sum += erlang.eval(state_age) * m_flows[i][idx_incoming_flow];
                }
                subcompartments[j] = 1. / (num_subcompartments * transition_rate) * sum;
                if (subcompartments[j] < 0) {
                    log_error(
                        "Initialization failed. Result for at least one subcompartment is less than zero. Please check "
                        "the input values.");
                    subcompartments[j] = -1.;
                    return subcompartments;
                }
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