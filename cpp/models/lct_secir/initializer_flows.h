/*
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/utils/metaprogramming.h"
#include "memilio/epidemiology/state_age_function.h"

namespace mio
{
namespace lsecir
{

/**
 * @brief Class that can be used to compute an initialization vector out of flows for an LCT Model
 *  with division in groups.
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
template <typename FP, typename Model>
class Initializer
{
public:
    using LctStatesGroups = typename Model::LctStatesGroups;

    /**
     * @brief Constructs a new Initializer object.
     *
     * @param[in] flows Initializing TimeSeries with flows fitting to these defined in InfectionTransition.
     *      For each group of m_model, InfectionTransition::Count entries are required.
     *      Timesteps should be equidistant and the values should be non-negative.
     *      The time history has to be long enough so that it is possible to calculate the initial vector.
     *      The length of the required time history depends on the Erlang densities used to compute the initial vector.
     * @param[in, out] model The LCT-SECIR model for which the initialization should be performed.
     */
    Initializer(TimeSeries<FP>&& flows, Model& model)
        : m_flows(std::move(flows))
        , m_model(model)
    {
        m_dt = m_flows.get_time(1) - m_flows.get_time(0);
    }

    /**
     * @brief Core function of Initializer.
     *
     * Computes a vector that can be used for the initialization of an LCT model stratified by groups with the number
     * of persons for each subcompartment for each group.
     * The initial value vector is updated in the model.
     *
     * @param[in] total_population The total size of the considered population.
     * @param[in] deaths Number of deceased people from the disease at time 0.
     * @param[in] total_confirmed_cases Total number of confirmed cases at time 0.
     * @return Returns true if one (or more) constraint(s) of the model, the initial flows
     *      or the computed initial value vector are not satisfied, otherwise false.
     */
    bool compute_initialization_vector(Eigen::VectorX<FP> const& total_population, Eigen::VectorX<FP> const& deaths,
                                       Eigen::VectorX<FP> const& total_confirmed_cases)
    {

        Eigen::VectorX<FP> init(m_model.populations.get_compartments());

        bool error = compute_initialization_vector_impl(init, total_population, deaths, total_confirmed_cases);
        if (error) {
            return true;
        }

        // Update initial value vector of the model.
        for (size_t i = 0; i < m_model.populations.get_num_compartments(); i++) {
            m_model.populations[i] = init[i];
        }

        // Check if constraints are fulfilled.
        return check_constraints();
    }

    /**
     * @brief Setter for the tolerance used to calculate the maximum support of ErlangDensity%s.
     *
     * @param[in] new_tol New tolerance.
     */
    void set_tol_for_support_max(FP new_tol)
    {
        m_tol = new_tol;
    }

private:
    /**
     * @brief Implementation of the calculation of the initial value vector slice that corresponds to a specified group.
     *
     * Computes a vector that can be used for the initialization of an LCT model stratified by groups with the number
     * of persons for each subcompartment. The groups are calculated recursively.
     *
     * @tparam Group The group for which the corresponding slice of the initial value vector is calculated.
     * @param[out] init The initial value vector under consideration.
     * @param[in] total_population The total size of the considered population.
     * @param[in] deaths Number of deceased people from the disease at time 0.
     * @param[in] total_confirmed_cases Total number of confirmed cases at time 0.
     * @return Returns true if one (or more) constraint(s) of the computed initial value vector are not satisfied,
     *       otherwise false.
     */
    template <size_t Group = 0>
    bool compute_initialization_vector_impl(Eigen::VectorX<FP>& init, Eigen::VectorX<FP> const& total_population,
                                            Eigen::VectorX<FP> const& deaths,
                                            Eigen::VectorX<FP> const& total_confirmed_cases)
    {
        static_assert((Group < Model::num_groups) && (Group >= 0), "The template parameter Group should be valid.");
        using LctStateGroup            = type_at_index_t<Group, LctStatesGroups>;
        Eigen::Index first_index       = m_model.populations.template get_first_index_of_group<Group>();
        Eigen::Index first_index_flows = Group * (Eigen::Index)InfectionTransition::Count;

        bool error =
            // Exposed.
            (compute_compartment<InfectionState::Exposed, Group>(
                init, first_index_flows + (Eigen::Index)InfectionTransition::SusceptibleToExposed,
                1 / m_model.parameters.template get<TimeExposed<FP>>()[Group])) ||
            // InfectedNoSymptoms.
            (compute_compartment<InfectionState::InfectedNoSymptoms, Group>(
                init, first_index_flows + (Eigen::Index)InfectionTransition::ExposedToInfectedNoSymptoms,
                1 / m_model.parameters.template get<TimeInfectedNoSymptoms<FP>>()[Group])) ||
            // InfectedSymptoms.
            (compute_compartment<InfectionState::InfectedSymptoms, Group>(
                init, first_index_flows + (Eigen::Index)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
                1 / m_model.parameters.template get<TimeInfectedSymptoms<FP>>()[Group])) ||
            // InfectedSevere.
            (compute_compartment<InfectionState::InfectedSevere, Group>(
                init, first_index_flows + (Eigen::Index)InfectionTransition::InfectedSymptomsToInfectedSevere,
                1 / m_model.parameters.template get<TimeInfectedSevere<FP>>()[Group])) ||
            // InfectedCritical.
            (compute_compartment<InfectionState::InfectedCritical, Group>(
                init, first_index_flows + (Eigen::Index)InfectionTransition::InfectedSevereToInfectedCritical,
                1 / m_model.parameters.template get<TimeInfectedCritical<FP>>()[Group]));
        if (error) {
            return true;
        }

        // Recovered.
        // Number of recovered is equal to the cumulative number of confirmed cases minus the number of people
        // in the group who are infected at the moment.
        init[first_index + LctStateGroup::template get_first_index<InfectionState::Recovered>()] =
            total_confirmed_cases[Group] -
            init.segment(first_index + LctStateGroup::template get_first_index<InfectionState::InfectedSymptoms>(),
                         LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms>() +
                             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere>() +
                             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>())
                .sum() -
            deaths[Group];

        // Susceptibles.
        init[first_index + LctStateGroup::template get_first_index<InfectionState::Susceptible>()] =
            total_population[Group] -
            init.segment(first_index + LctStateGroup::template get_first_index<InfectionState::Exposed>(),
                         LctStateGroup::Count - 2)
                .sum() -
            deaths[Group];

        // Dead.
        init[first_index + LctStateGroup::template get_first_index<InfectionState::Dead>()] = deaths[Group];

        for (size_t state_idx : {LctStateGroup::template get_first_index<InfectionState::Susceptible>(),
                                 LctStateGroup::template get_first_index<InfectionState::Recovered>(),
                                 LctStateGroup::template get_first_index<InfectionState::Dead>()}) {
            if (floating_point_less<FP>(init[first_index + state_idx], 0.0, 1e-8)) {
                log_error("Initialization failed. Values of total_confirmed_cases, deaths and total_population do not "
                          "match the specified values for the transitions, leading to at least one negative result.");
                return true;
            }
        }
        // Function call for next group if applicable.
        if constexpr (Group + 1 < Model::num_groups) {
            return compute_initialization_vector_impl<Group + 1>(init, total_population, deaths, total_confirmed_cases);
        }
        else {
            return false;
        }
    }

    /**
     * @brief Checks constraints of the Initializer including checks for the model.
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false.
     */
    bool check_constraints() const
    {
        if (!((Eigen::Index)InfectionTransition::Count * Model::num_groups == m_flows.get_num_elements())) {
            log_error("Initial condition size does not match subcompartments.");
            return true;
        }

        if (!(floating_point_equal<FP>(0., m_flows.get_last_time(), 1e-8, 1e-14))) {
            log_error("Last time point in flows has to be 0.");
            return true;
        }

        for (int i = 1; i < m_flows.get_num_time_points(); i++) {
            if (!(floating_point_equal<FP>(m_dt, m_flows.get_time(i) - m_flows.get_time(i - 1), 1e-8, 1e-14))) {
                log_error("Time points in flows have to be equidistant.");
                return true;
            }
        }

        if (!(m_dt < 1)) {
            log_warning("Step size was set very large. The result could be distorted.");
            return true;
        }

        // Check if model is valid and the calculated initial value vector is valid.
        return m_model.check_constraints();
    }

    /**
     * @brief Computes a slice of the initial value vector for each subcompartment of one InfectionState
     *  for a specified group.
     *
     * @tparam Group The group for which the corresponding slice of the initial value vector is calculated.
     * @tparam State The InfectionState for which the corresponding slice of the initial value vector is calculated
     * @param[out] init The initial value vector under consideration.
     * @param[in] idx_incoming_flow Index of the flow which is relevant for the calculation, so the flow
     *      to the InfectionState.
     * @param[in] transition_rate Specifies the transition rate of the InfectionState.
     *      'This is equal to 1 / (expected time in the InfectionState).

     * @return Returns true if one (or more) constraint(s) of the computed slice of the initial value vector
     *      are not satisfied, otherwise false.
     */
    template <InfectionState State, size_t Group>
    bool compute_compartment(Eigen::VectorX<FP>& init, Eigen::Index idx_incoming_flow, FP transition_rate) const
    {
        using std::ceil;

        size_t first_index = m_model.populations.template get_first_index_of_group<Group>() +
                             type_at_index_t<Group, LctStatesGroups>::template get_first_index<State>();
        size_t num_subcompartments = type_at_index_t<Group, LctStatesGroups>::template get_num_subcompartments<State>();

        // Initialize relevant density for the InfectionState.
        // For the first subcompartment a shape parameter of one is needed.
        ErlangDensity<FP> erlang(1, 1. / (num_subcompartments * transition_rate));

        // Initialize other relevant parameters.
        FP calc_time{0};
        Eigen::Index calc_time_index{0};
        FP state_age{0};
        FP sum{0};
        Eigen::Index num_time_points = m_flows.get_num_time_points();

        // Calculate number of people in each subcompartment.
        for (size_t j = 0; j < num_subcompartments; j++) {
            // For subcompartment number j+1, shape parameter j+1 is needed.
            erlang.set_distribution_parameter((FP)j + 1);

            // Determine relevant calculation area and corresponding index.
            calc_time       = erlang.get_support_max(m_dt, m_tol);
            calc_time_index = (Eigen::Index)ceil(calc_time / m_dt) - 1;

            if (num_time_points < calc_time_index) {
                log_error("Initialization failed. Not enough time points for the transitions are given. More than {} "
                          "are needed but just {} are given.",
                          calc_time_index, num_time_points);
                return true;
            }
            else {
                // Approximate integral with non-standard scheme.
                for (Eigen::Index i = num_time_points - calc_time_index; i < num_time_points; i++) {
                    state_age = (num_time_points - i) * m_dt;
                    sum += erlang.eval(state_age) * m_flows[i][idx_incoming_flow];
                }
                init[first_index + j] = 1. / (num_subcompartments * transition_rate) * sum;
                if (init[first_index + j] < 0) {
                    log_error(
                        "Initialization failed. Result for at least one subcompartment is less than zero. Please check "
                        "the input values.");
                    return true;
                }
            }
            sum = 0;
        }
        return false;
    }

    TimeSeries<FP> m_flows; ///< TimeSeries with the flows which are used to calculate the initial vector.
    Model& m_model; ///< The LCT-SECIR model for which the initialization should be performed.
    FP m_dt{}; ///< Step size of the times in m_flows and time step for the approximation of the integral.
    FP m_tol{1e-10}; ///< Tolerance used to calculate the maximum support of the ErlangDensity%s.
};

} // namespace lsecir
} // namespace mio

#endif
