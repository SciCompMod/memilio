/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Annika Jungklaus, Lena Ploetzke
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

#ifndef LCT_SECIR_2_DISEASE_MODEL_H
#define LCT_SECIR_2_DISEASE_MODEL_H

#include "lct_secir_2_diseases/parameters.h"
#include "lct_secir_2_diseases/infection_state.h"
#include "memilio/compartments/compartmental_model.h"
#include "memilio/epidemiology/lct_populations.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/type_list.h"
#include "memilio/utils/metaprogramming.h"

namespace mio
{
namespace lsecir2d
{
/**
 * @brief Class that defines an LCT-SECIR-2-DISEASE model.
 *
 * @tparam LctStates The LCT2D model can work with any number of LctStates, where each LctState corresponds to a group,
 * e.g. one AgeGroup. The purpose of the LctStates is to define the number of subcompartments for each InfectionState.
 * If you do not want to divide the population into groups, just use one LctState.
 * If you want to divide the population according to more than one category, e.g. sex and age, 
 * you have to specify one LctState for each pair of groups, e.g. for (female, A00-A04), (female, A05-A14) etc. 
 * This is because the number of subcompartments can be different for each group.
 * Therefore, the number of LctStates also determines the number of groups. 
 */
template <typename FP, class... LctStates>
class Model : public CompartmentalModel<FP, InfectionState, LctPopulations<FP, LctStates...>, Parameters<FP>>
{
public:
    using LctStatesGroups = TypeList<LctStates...>;
    using Base            = CompartmentalModel<FP, InfectionState, LctPopulations<FP, LctStates...>, Parameters<FP>>;
    using typename Base::ParameterSet;
    using typename Base::Populations;
    static size_t constexpr num_groups = sizeof...(LctStates);
    static_assert(num_groups >= 1, "The number of LctStates provided should be at least one.");

    /// @brief Default constructor.
    Model()
        : Base(Populations(), ParameterSet(num_groups))
    {
    }

    /** @brief Constructor using Populations and ParameterSet.
     * @param[in] pop An instance of the Populations class.
     * @param[in] params Parameters used to be used in the simulation.
     */
    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    /**
     * @brief Evaluates the right-hand-side f of the ODE dydt = f(y, t).
     *
     * The LCT-SECIR-2-DISEASE model is defined through ordinary differential equations of the form dydt = f(y, t). 
     * y is a vector containing the number of individuals for each (sub-) compartment.
     * This function evaluates the right-hand-side f of the ODE and can be used in an ODE solver.
     * @param[in] pop The current state of the population in the geographic unit we are considering.
     * @param[in] y The current state of the model (or a subpopulation) as a flat array.
     * @param[in] t The current time.
     * @param[out] dydt A reference to the calculated output.
     */
    void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                         Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
    {
        // Vectors are sorted such that we first have all InfectionState%s for AgeGroup 0,
        // afterwards all for AgeGroup 1 and so on.
        dydt.setZero();
        get_derivatives_impl(pop, y, t, dydt);
    }

    /**
     * @brief Cumulates a simulation result with subcompartments to produce a result that divides the population only
     *   into the infection states defined in InfectionState.
     *
     * If the model is used for simulation, we will get a result in form of a TimeSeries with infection states divided 
     * in subcompartments.
     * The function calculates a TimeSeries without subcompartments from another TimeSeries with subcompartments. 
     * This is done by summing up the numbers in the subcompartments.
     * @param[in] subcompartments_ts Result of a simulation with the model.
     * @return Result of the simulation divided in infection states without subcompartments. 
     *  Returns TimeSeries with values -1 if calculation is not possible.
     */
    TimeSeries<FP> calculate_compartments(const TimeSeries<FP>& subcompartments_ts) const
    {
        Eigen::Index count_InfStates  = (Eigen::Index)InfectionState::Count;
        Eigen::Index num_compartments = count_InfStates * num_groups;
        TimeSeries<FP> compartments_ts(num_compartments);
        if (!(this->populations.get_num_compartments() == (size_t)subcompartments_ts.get_num_elements())) {
            log_error("Result does not match InfectionState of the Model.");
            Eigen::VectorX<FP> wrong_size = Eigen::VectorX<FP>::Constant(num_compartments, -1);
            compartments_ts.add_time_point(-1, wrong_size);
            return compartments_ts;
        }
        Eigen::VectorX<FP> compartments(num_compartments);
        for (Eigen::Index timepoint = 0; timepoint < subcompartments_ts.get_num_time_points(); ++timepoint) {
            compress_vector(subcompartments_ts[timepoint], compartments);
            compartments_ts.add_time_point(subcompartments_ts.get_time(timepoint), compartments);
        }

        return compartments_ts;
    }

    /**
     * @brief Checks that the model satisfies all constraints (e.g. parameter or population constraints).
     * @return Returns true if one or more constraints are not satisfied, false otherwise.
     */
    bool check_constraints() const
    {

        return (check_constraints_impl() || this->parameters.check_constraints() ||
                this->populations.check_constraints());
    }

private:
    /**
     * @brief Converts a vector with subcompartments in a vector without subcompartments 
     * by summing up subcompartment values.
     * This is done recursively for each group which corresponds to a slice of the vector.
     *
     * @tparam group The group specifying the slice of the vector being considered. 
     * @param[in] subcompartments The vector that should be converted.
     * @param[out] compartments Reference to the vector where the output is stored.
     */
    template <size_t Group = 0>
    void compress_vector(const Eigen::VectorX<FP>& subcompartments, Eigen::VectorX<FP>& compartments) const
    {
        static_assert((Group < num_groups) && (Group >= 0), "The template parameter Group should be valid.");
        using LctStateGroup = type_at_index_t<Group, LctStates...>;

        // Define first index of the group Group in a vector including all compartments without a resolution
        // in subcompartments.
        Eigen::Index count_InfStates         = (Eigen::Index)InfectionState::Count;
        Eigen::Index first_index_group_comps = Group * count_InfStates;

        // Use function from the LctState of the Group to calculate the vector without subcompartments
        // using the corresponding vector with subcompartments.
        compartments.segment(first_index_group_comps, count_InfStates) =
            LctStateGroup::calculate_compartments(subcompartments.segment(
                this->populations.template get_first_index_of_group<Group>(), LctStateGroup::Count));

        // Function call for next group if applicable.
        if constexpr (Group + 1 < num_groups) {
            compress_vector<Group + 1>(subcompartments, compartments);
        }
    }

    /**
     * @brief Evaluates the right-hand-side f of the ODE dydt = f(y, t) recursively for each group.
     *
     * See also the function get_derivative.
     * For each group, one slice of the output vector is calculated.
     * @tparam group The group specifying the slice of the vector being considered.  
     * @param[in] pop The current state of the population in the geographic unit we are considering.
     * @param[in] y The current state of the model (or a subpopulation) as a flat array.
     * @param[in] t The current time.
     * @param[out] dydt A reference to the calculated output.
     */
    template <size_t Group = 0>
    void get_derivatives_impl(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                              Eigen::Ref<Eigen::VectorX<FP>> dydt) const
    {
        static_assert((Group < num_groups) && (Group >= 0), "The template parameter Group should be valid.");
        using LctStateGroup = type_at_index_t<Group, LctStates...>;

        size_t first_index_group = this->populations.template get_first_index_of_group<Group>();
        const auto& params       = this->parameters;
        FP flow                  = 0;

        // Indices of first subcompartment of the InfectionState for the group in the vectors.
        size_t Ei_1a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::Exposed_1a>();
        size_t Ei_2a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::Exposed_2a>();
        size_t Ei_1b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::Exposed_1b>();
        size_t Ei_2b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::Exposed_2b>();
        size_t INSi_1a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedNoSymptoms_1a>();
        size_t INSi_2a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedNoSymptoms_2a>();
        size_t INSi_1b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedNoSymptoms_1b>();
        size_t INSi_2b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedNoSymptoms_2b>();
        size_t ISyi_1a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSymptoms_1a>();
        size_t ISyi_2a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSymptoms_2a>();
        size_t ISyi_1b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSymptoms_1b>();
        size_t ISyi_2b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSymptoms_2b>();
        size_t ISevi_1a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSevere_1a>();
        size_t ISevi_2a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSevere_2a>();
        size_t ISevi_1b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSevere_1b>();
        size_t ISevi_2b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSevere_2b>();
        size_t ICri_1a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedCritical_1a>();
        size_t ICri_2a_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedCritical_2a>();
        size_t ICri_1b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedCritical_1b>();
        size_t ICri_2b_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedCritical_2b>();
        size_t Ri_a  = first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered_1a>();
        size_t Ri_b  = first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered_1b>();
        size_t Ri_ab = first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered_2ab>();
        size_t Di_a  = first_index_group + LctStateGroup::template get_first_index<InfectionState::Dead_a>();
        size_t Di_b  = first_index_group + LctStateGroup::template get_first_index<InfectionState::Dead_b>();

        // Calculate derivative of the Susceptible compartment.
        // The outflow is generated by disease a and disease b both.
        FP part_a = 0.; // Part of people getting infected with disease a.
        FP part_b = 0.; // Part of people getting infected with disease b.
        interact<Group, 0>(pop, y, t, dydt, &part_a, &part_b, 2);
        // Split flow.
        FP div_part_both = ((part_a + part_b) < Limits<FP>::zero_tolerance()) ? 0.0 : 1.0 / (part_a + part_b);

        // Start with derivatives of first infections (X_1a, X_1b)
        // Calculate derivative of the Exposed_1x compartments, split flow from Susceptible.
        // Exposed 1 a:
        dydt[Ei_1a_first_index] = -dydt[first_index_group] * part_a * div_part_both;
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_1a>(); subcomp++) {
            // Variable flow stores the value of the flow from one subcompartment to the next one.
            // Ei_1a_first_index + subcomp is always the index of a (sub-)compartment of Exposed_1a and Ei_1a_first_index
            // + subcomp + 1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_1a>() *
                   (1 / params.template get<TimeExposed_a<FP>>()[Group]) * y[Ei_1a_first_index + subcomp];
            // Subtract flow from dydt[Ei_1a_first_index + subcomp] and add to next subcompartment.
            dydt[Ei_1a_first_index + subcomp] -= flow;
            dydt[Ei_1a_first_index + subcomp + 1] = flow;
        }
        // Exposed 1 b:
        dydt[Ei_1b_first_index] = -dydt[first_index_group] * part_b * div_part_both;
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_1b>(); subcomp++) {
            // Variable flow stores the value of the flow from one subcompartment to the next one.
            // Ei_1a_first_index + subcomp is always the index of a (sub-)compartment of Exposed_1a and Ei_1b_first_index
            // + subcomp + 1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_1b>() *
                   (1 / params.template get<TimeExposed_b<FP>>()[Group]) * y[Ei_1b_first_index + subcomp];
            // Subtract flow from dydt[Ei_1b_first_index + subcomp] and add to next subcompartment.
            dydt[Ei_1b_first_index + subcomp] -= flow;
            dydt[Ei_1b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedNoSymptoms_1a compartment
        // flow from each subcompartment to the next.
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1a>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1a>() *
                   (1 / params.template get<TimeInfectedNoSymptoms_a<FP>>()[Group]) * y[INSi_1a_first_index + subcomp];
            dydt[INSi_1a_first_index + subcomp] -= flow;
            dydt[INSi_1a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedNoSymptoms_1b compartment
        // flow from each subcompartment to the next.
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1b>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1b>() *
                   (1 / params.template get<TimeInfectedNoSymptoms_b<FP>>()[Group]) * y[INSi_1b_first_index + subcomp];
            dydt[INSi_1b_first_index + subcomp] -= flow;
            dydt[INSi_1b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSymptoms_1a compartment.
        // Flow from last (sub-) compartment of InfectedNoSymptoms_1a must be split between
        // the first subcompartment of InfectedSymptoms_1a and Recovered_1a.
        dydt[Ri_a] = dydt[ISyi_1a_first_index] * params.template get<RecoveredPerInfectedNoSymptoms_a<FP>>()[Group];
        dydt[ISyi_1a_first_index] =
            dydt[ISyi_1a_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms_a<FP>>()[Group]);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_1a>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_1a>() *
                   (1 / params.template get<TimeInfectedSymptoms_a<FP>>()[Group]) * y[ISyi_1a_first_index + subcomp];
            dydt[ISyi_1a_first_index + subcomp] -= flow;
            dydt[ISyi_1a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedSymptoms_1b compartment.
        // Flow from last (sub-) compartment of InfectedNoSymptoms_1b must be split between
        // the first subcompartment of InfectedSymptoms_1b and Recovered_1b.
        dydt[Ri_b] = dydt[ISyi_1b_first_index] * params.template get<RecoveredPerInfectedNoSymptoms_b<FP>>()[Group];
        dydt[ISyi_1b_first_index] =
            dydt[ISyi_1b_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms_b<FP>>()[Group]);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_1b>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_1b>() *
                   (1 / params.template get<TimeInfectedSymptoms_b<FP>>()[Group]) * y[ISyi_1b_first_index + subcomp];
            dydt[ISyi_1b_first_index + subcomp] -= flow;
            dydt[ISyi_1b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSevere_1a compartment.
        // Again split the flow from the last subcompartment of InfectedSymptoms_1a.
        dydt[Ri_a] += dydt[ISevi_1a_first_index] * (1 - params.template get<SeverePerInfectedSymptoms_a<FP>>()[Group]);
        dydt[ISevi_1a_first_index] =
            dydt[ISevi_1a_first_index] * params.template get<SeverePerInfectedSymptoms_a<FP>>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_1a>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_1a>() *
                   (1 / params.template get<TimeInfectedSevere_a<FP>>()[Group]) * y[ISevi_1a_first_index + subcomp];
            dydt[ISevi_1a_first_index + subcomp] -= flow;
            dydt[ISevi_1a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedSevere_1b compartment.
        // Again split the flow from the last subcompartment of InfectedSymptoms_1b.
        dydt[Ri_b] += dydt[ISevi_1b_first_index] * (1 - params.template get<SeverePerInfectedSymptoms_b<FP>>()[Group]);
        dydt[ISevi_1b_first_index] =
            dydt[ISevi_1b_first_index] * params.template get<SeverePerInfectedSymptoms_b<FP>>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_1b>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_1b>() *
                   (1 / params.template get<TimeInfectedSevere_b<FP>>()[Group]) * y[ISevi_1b_first_index + subcomp];
            dydt[ISevi_1b_first_index + subcomp] -= flow;
            dydt[ISevi_1b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedCritical_1a compartment.
        // Again split flow from last subcompartment of InfectedSevere_1a between Recovered_1a and InfectedCritical_1a.
        dydt[Ri_a] += dydt[ICri_1a_first_index] * (1 - params.template get<CriticalPerSevere_a<FP>>()[Group]);
        dydt[ICri_1a_first_index] = dydt[ICri_1a_first_index] * params.template get<CriticalPerSevere_a<FP>>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() - 1;
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() *
                   (1 / params.template get<TimeInfectedCritical_a<FP>>()[Group]) * y[ICri_1a_first_index + subcomp];
            dydt[ICri_1a_first_index + subcomp] -= flow;
            dydt[ICri_1a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedCritical_1b compartment.
        // Again split flow from last subcompartment of InfectedSevere_1b between Recovered_1b and InfectedCritical_1b.
        dydt[Ri_b] += dydt[ICri_1b_first_index] * (1 - params.template get<CriticalPerSevere_b<FP>>()[Group]);
        dydt[ICri_1b_first_index] = dydt[ICri_1b_first_index] * params.template get<CriticalPerSevere_b<FP>>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() - 1;
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() *
                   (1 / params.template get<TimeInfectedCritical_b<FP>>()[Group]) * y[ICri_1b_first_index + subcomp];
            dydt[ICri_1b_first_index + subcomp] -= flow;
            dydt[ICri_1b_first_index + subcomp + 1] = flow;
        }

        // Flow from InfectedCritical compartments has to be divided between Recovered and Dead compartments.
        // Must be calculated separately in order not to overwrite the already calculated values ​​for Recovered.
        // Outflow from InfectedCritical_1a:
        flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() *
               (1 / params.template get<TimeInfectedCritical_a<FP>>()[Group]) *
               y[ICri_1a_first_index +
                 LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() - 1];
        dydt[ICri_1a_first_index +
             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() - 1] -= flow;
        dydt[Ri_a] = dydt[Ri_a] + (1 - params.template get<DeathsPerCritical_a<FP>>()[Group]) * flow;
        dydt[Di_a] = dydt[Di_a] + params.template get<DeathsPerCritical_a<FP>>()[Group] * flow;
        // Outflow from InfectedCritical_1b:
        flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() *
               (1 / params.template get<TimeInfectedCritical_b<FP>>()[Group]) *
               y[ICri_1b_first_index +
                 LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() - 1];
        dydt[ICri_1b_first_index +
             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() - 1] -= flow;
        dydt[Ri_b] = dydt[Ri_b] + (1 - params.template get<DeathsPerCritical_b<FP>>()[Group]) * flow;
        dydt[Di_b] = dydt[Di_b] + params.template get<DeathsPerCritical_b<FP>>()[Group] * flow;

        // Second Infection:
        // Outflow from Recovered_1a and Recovered_1b (people getting infected for the second time).
        double temp_Ra = dydt[Ri_a];
        // Outflow from Recovered_1a is only affected by b.
        interact<Group, 0>(pop, y, t, dydt, &part_a, &part_b, 1);
        double temp_Rb = dydt[Ri_b];
        // Outflow from Recovered_1b is only affected by a.
        interact<Group, 0>(pop, y, t, dydt, &part_a, &part_b, 0);

        // Calculate derivative of the Exposed_2i compartments:
        // Exposed_2a:
        dydt[Ei_2a_first_index] = -(dydt[Ri_b] - temp_Rb);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_2a>(); subcomp++) {
            // Variable flow stores the value of the flow from one subcompartment to the next one.
            // Ei_2a_first_index + subcomp is always the index of a (sub-)compartment of Exposed and Ei_2a_first_index
            // + subcomp + 1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_2a>() *
                   (1 / params.template get<TimeExposed_a<FP>>()[Group]) * y[Ei_2a_first_index + subcomp];
            // Subtract flow from dydt[Ei_2a_first_index + subcomp] and add to next subcompartment.
            dydt[Ei_2a_first_index + subcomp] -= flow;
            dydt[Ei_2a_first_index + subcomp + 1] = flow;
        }
        // Exposed_2b:
        dydt[Ei_2b_first_index] = -(dydt[Ri_a] - temp_Ra);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_2b>(); subcomp++) {
            // Variable flow stores the value of the flow from one subcompartment to the next one.
            // Ei_2b_first_index + subcomp is always the index of a (sub-)compartment of Exposed and Ei_2b_first_index
            // + subcomp + 1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_2b>() *
                   (1 / params.template get<TimeExposed_b<FP>>()[Group]) * y[Ei_2b_first_index + subcomp];
            // Subtract flow from dydt[Ei_2b_first_index + subcomp] and add to next subcompartment.
            dydt[Ei_2b_first_index + subcomp] -= flow;
            dydt[Ei_2b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedNoSymptoms_2a compartment.
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2a>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2a>() *
                   (1 / params.template get<TimeInfectedNoSymptoms_a<FP>>()[Group]) * y[INSi_2a_first_index + subcomp];
            dydt[INSi_2a_first_index + subcomp] -= flow;
            dydt[INSi_2a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedNoSymptoms_2b compartment.
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2b>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2b>() *
                   (1 / params.template get<TimeInfectedNoSymptoms_b<FP>>()[Group]) * y[INSi_2b_first_index + subcomp];
            dydt[INSi_2b_first_index + subcomp] -= flow;
            dydt[INSi_2b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSymptoms_2a compartment.
        // Flow from last (sub-) compartment of InfectedNoSymptoms_2a must be split between
        // the first subcompartment of InfectedSymptoms_2a and Recovered_2ab.
        dydt[Ri_ab] += dydt[ISyi_2a_first_index] * params.template get<RecoveredPerInfectedNoSymptoms_a<FP>>()[Group];
        dydt[ISyi_2a_first_index] =
            dydt[ISyi_2a_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms_a<FP>>()[Group]);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_2a>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_2a>() *
                   (1 / params.template get<TimeInfectedSymptoms_a<FP>>()[Group]) * y[ISyi_2a_first_index + subcomp];
            dydt[ISyi_2a_first_index + subcomp] -= flow;
            dydt[ISyi_2a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedSymptoms_2b compartment.
        // Flow from last (sub-) compartment of InfectedNoSymptoms_2b must be split between
        // the first subcompartment of InfectedSymptoms_2b and Recovered_2ab.
        dydt[Ri_ab] += dydt[ISyi_2b_first_index] * params.template get<RecoveredPerInfectedNoSymptoms_b<FP>>()[Group];
        dydt[ISyi_2b_first_index] =
            dydt[ISyi_2b_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms_b<FP>>()[Group]);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_2b>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_2b>() *
                   (1 / params.template get<TimeInfectedSymptoms_b<FP>>()[Group]) * y[ISyi_2b_first_index + subcomp];
            dydt[ISyi_2b_first_index + subcomp] -= flow;
            dydt[ISyi_2b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSevere_2a compartment.
        // Again split the flow from the last subcompartment of InfectedSymptoms_2a.
        dydt[Ri_ab] += dydt[ISevi_2a_first_index] * (1 - params.template get<SeverePerInfectedSymptoms_a<FP>>()[Group]);
        dydt[ISevi_2a_first_index] =
            dydt[ISevi_2a_first_index] * params.template get<SeverePerInfectedSymptoms_a<FP>>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_2a>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_2a>() *
                   (1 / params.template get<TimeInfectedSevere_a<FP>>()[Group]) * y[ISevi_2a_first_index + subcomp];
            dydt[ISevi_2a_first_index + subcomp] -= flow;
            dydt[ISevi_2a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedSevere compartment.
        // Again split the flow from the last subcompartment of InfectedSymptoms_2b.
        dydt[Ri_ab] += dydt[ISevi_2b_first_index] * (1 - params.template get<SeverePerInfectedSymptoms_b<FP>>()[Group]);
        dydt[ISevi_2b_first_index] =
            dydt[ISevi_2b_first_index] * params.template get<SeverePerInfectedSymptoms_b<FP>>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_2b>();
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_2b>() *
                   (1 / params.template get<TimeInfectedSevere_b<FP>>()[Group]) * y[ISevi_2b_first_index + subcomp];
            dydt[ISevi_2b_first_index + subcomp] -= flow;
            dydt[ISevi_2b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedCritical compartment.
        // Again split flow from last subcompartment of InfectedSevere_2a between Recovered_2ab and InfectedCritical_2a.
        dydt[Ri_ab] += dydt[ICri_2a_first_index] * (1 - params.template get<CriticalPerSevere_a<FP>>()[Group]);
        dydt[ICri_2a_first_index] = dydt[ICri_2a_first_index] * params.template get<CriticalPerSevere_a<FP>>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() - 1;
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() *
                   (1 / params.template get<TimeInfectedCritical_a<FP>>()[Group]) * y[ICri_2a_first_index + subcomp];
            dydt[ICri_2a_first_index + subcomp] -= flow;
            dydt[ICri_2a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedCritical compartment.
        // Again split flow from last subcompartment of InfectedSevere_2b between Recovered_2ab and InfectedCritical_2b.
        dydt[Ri_ab] += dydt[ICri_2b_first_index] * (1 - params.template get<CriticalPerSevere_b<FP>>()[Group]);
        dydt[ICri_2b_first_index] = dydt[ICri_2b_first_index] * params.template get<CriticalPerSevere_b<FP>>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() - 1;
             subcomp++) {
            flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() *
                   (1 / params.template get<TimeInfectedCritical_b<FP>>()[Group]) * y[ICri_2b_first_index + subcomp];
            dydt[ICri_2b_first_index + subcomp] -= flow;
            dydt[ICri_2b_first_index + subcomp + 1] = flow;
        }

        // Last flow from InfectedCritical has to be divided between Recovered and Dead.
        // Must be calculated separately in order not to overwrite the already calculated values ​​for Recovered.
        // Outflow from InfectedCritical_2a:
        flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() *
               (1 / params.template get<TimeInfectedCritical_a<FP>>()[Group]) *
               y[ICri_2a_first_index +
                 LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() - 1];
        dydt[ICri_2a_first_index +
             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() - 1] -= flow;
        dydt[Ri_ab] = dydt[Ri_ab] + (1 - params.template get<DeathsPerCritical_a<FP>>()[Group]) * flow;
        dydt[Di_a]  = dydt[Di_a] + params.template get<DeathsPerCritical_a<FP>>()[Group] * flow;
        // Outflow from InfectedCritical_2b:
        flow = (FP)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() *
               (1 / params.template get<TimeInfectedCritical_b<FP>>()[Group]) *
               y[ICri_2b_first_index +
                 LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() - 1];
        dydt[ICri_2b_first_index +
             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() - 1] -= flow;
        dydt[Ri_ab] = dydt[Ri_ab] + (1 - params.template get<DeathsPerCritical_b<FP>>()[Group]) * flow;
        dydt[Di_b]  = dydt[Di_b] + params.template get<DeathsPerCritical_b<FP>>()[Group] * flow;

        // Function call for next group if applicable.
        if constexpr (Group + 1 < num_groups) {
            get_derivatives_impl<Group + 1>(pop, y, t, dydt);
        }
    }

    /**
     * @brief Calculates flows that are caused by people becoming infected (outflow from compartment S, Ra or Rb) for Group1.
     *
     * This is done recursively by calculating the interaction terms with each group.
     * @tparam Group1 The group for which the derivative of the compartment should be calculated. 
     * @tparam Group2 The group that Group1 interacts with. 
     * @param[in] pop The current state of the population in the geographic unit we are considering.
     * @param[in] y The current state of the model (or a subpopulation) as a flat array.
     * @param[in] t The current time.
     * @param[out] dydt A reference to the calculated output.
     * The flow from Susceptible needs to be split into 2 parts for Exposed_1a and Exposed_1b:
     * @param[out] part_a Reference to amount of flow caused by people infected with disease a.
     * @param[out] part_b Reference to amount of flow caused by people infected with disease b.
     * @param[in] relevant_disease Index for which infected people cause the outflow (0 = a, 1 = b, 2 = a and b).
     */
    template <size_t Group1, size_t Group2 = 0>
    void interact(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                  Eigen::Ref<Eigen::VectorX<FP>> dydt, double* part_a, double* part_b, int relevant_disease) const
    {
        static_assert((Group1 < num_groups) && (Group1 >= 0) && (Group2 < num_groups) && (Group2 >= 0),
                      "The template parameters Group1 & Group2 should be valid.");
        using LctStateGroup1           = type_at_index_t<Group1, LctStates...>;
        using LctStateGroup2           = type_at_index_t<Group2, LctStates...>;
        FP InfectedNoSymptoms_Group2_a = 0;
        FP InfectedSymptoms_Group2_a   = 0;
        FP InfectedNoSymptoms_Group2_b = 0;
        FP InfectedSymptoms_Group2_b   = 0;
        const auto& params             = this->parameters;

        size_t first_index_group1 = this->populations.template get_first_index_of_group<Group1>();
        size_t first_index_group2 = this->populations.template get_first_index_of_group<Group2>();

        // Calculate sum of all subcompartments for InfectedNoSymptoms for disease a of Group2.
        InfectedNoSymptoms_Group2_a =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedNoSymptoms_1a>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1a>())
                .sum() +
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedNoSymptoms_2a>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2a>())
                .sum();
        // Calculate sum of all subcompartments for InfectedSymptoms for disease a of Group2.
        InfectedSymptoms_Group2_a =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedSymptoms_1a>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedSymptoms_1a>())
                .sum() +
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedSymptoms_2a>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedSymptoms_2a>())
                .sum();
        // Calculate sum of all subcompartments for InfectedNoSymptoms for disease b of Group2.
        InfectedNoSymptoms_Group2_b =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedNoSymptoms_1b>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1b>())
                .sum() +
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedNoSymptoms_2b>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2b>())
                .sum();
        // Calculate sum of all subcompartments for InfectedSymptoms for disease b of Group2.
        InfectedSymptoms_Group2_b =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedSymptoms_1b>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedSymptoms_1b>())
                .sum() +
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedSymptoms_2b>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedSymptoms_2b>())
                .sum();
        // Size of the subpopulation Group2 without dead people.
        FP N_Group2 = pop.segment(first_index_group2, LctStateGroup2::Count).sum(); // Sum over all compartments.
        N_Group2 = N_Group2 - pop.segment(LctStateGroup2::template get_first_index<InfectionState::Dead_a>(), 1).sum() -
                   pop.segment(LctStateGroup2::template get_first_index<InfectionState::Dead_b>(), 1)
                       .sum(); // All people minus dead people.
        const FP div_N_Group2 = (N_Group2 < Limits<FP>::zero_tolerance()) ? 0.0 : 1.0 / N_Group2;
        FP season_val         = 1 + params.template get<Seasonality<FP>>() *
                                sin(3.141592653589793 * ((params.template get<StartDay<FP>>() + t) / 182.5 + 0.5));

        if (relevant_disease == 0) { // Disease a.
            // Get index for compartment Recovered_1b and calculate outflow driven by disease a.
            size_t compartment_index =
                first_index_group1 + LctStateGroup1::template get_first_index<InfectionState::Recovered_1b>();
            dydt[compartment_index] +=
                -y[compartment_index] * div_N_Group2 * season_val *
                params.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(SimulationTime<FP>(t))(
                    static_cast<Eigen::Index>(Group1), static_cast<Eigen::Index>(Group2)) *
                params.template get<TransmissionProbabilityOnContact_a<FP>>()[Group1] *
                (params.template get<RelativeTransmissionNoSymptoms_a<FP>>()[Group2] * InfectedNoSymptoms_Group2_a +
                 params.template get<RiskOfInfectionFromSymptomatic_a<FP>>()[Group2] * InfectedSymptoms_Group2_a);
        }
        else if (relevant_disease == 1) { // Disease b.
            // Get index for compartment Recovered_1a and calculate outflow driven by disease b.
            size_t compartment_index =
                first_index_group1 + LctStateGroup1::template get_first_index<InfectionState::Recovered_1a>();
            dydt[compartment_index] +=
                -y[compartment_index] * div_N_Group2 * season_val *
                params.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(SimulationTime<FP>(t))(
                    static_cast<Eigen::Index>(Group1), static_cast<Eigen::Index>(Group2)) *
                params.template get<TransmissionProbabilityOnContact_b<FP>>()[Group1] *
                (params.template get<RelativeTransmissionNoSymptoms_b<FP>>()[Group2] * InfectedNoSymptoms_Group2_b +
                 params.template get<RiskOfInfectionFromSymptomatic_b<FP>>()[Group2] * InfectedSymptoms_Group2_b);
        }
        else if (relevant_disease == 2) { // Both diseases drive outflow from the Susceptible compartment.
            size_t compartment_index = first_index_group1;
            dydt[compartment_index] +=
                -y[compartment_index] * div_N_Group2 * season_val *
                params.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(SimulationTime<FP>(t))(
                    static_cast<Eigen::Index>(Group1), static_cast<Eigen::Index>(Group2)) *
                (params.template get<TransmissionProbabilityOnContact_a<FP>>()[Group1] *
                     (params.template get<RelativeTransmissionNoSymptoms_a<FP>>()[Group2] *
                          InfectedNoSymptoms_Group2_a +
                      params.template get<RiskOfInfectionFromSymptomatic_a<FP>>()[Group2] * InfectedSymptoms_Group2_a) +
                 params.template get<TransmissionProbabilityOnContact_b<FP>>()[Group1] *
                     (params.template get<RelativeTransmissionNoSymptoms_b<FP>>()[Group2] *
                          InfectedNoSymptoms_Group2_b +
                      params.template get<RiskOfInfectionFromSymptomatic_b<FP>>()[Group2] * InfectedSymptoms_Group2_b));
        }
        // To split the outflow from Susceptible between Exposed_1a and Exposed_1b.
        *part_a += params.template get<TransmissionProbabilityOnContact_a<FP>>()[Group1] *
                   (params.template get<RelativeTransmissionNoSymptoms_a<FP>>()[Group2] * InfectedNoSymptoms_Group2_a +
                    params.template get<RiskOfInfectionFromSymptomatic_a<FP>>()[Group2] * InfectedSymptoms_Group2_a);
        *part_b += params.template get<TransmissionProbabilityOnContact_b<FP>>()[Group1] *
                   (params.template get<RelativeTransmissionNoSymptoms_b<FP>>()[Group2] * InfectedNoSymptoms_Group2_b +
                    params.template get<RiskOfInfectionFromSymptomatic_b<FP>>()[Group2] * InfectedSymptoms_Group2_b);

        if constexpr (Group2 + 1 < num_groups) {
            interact<Group1, Group2 + 1>(pop, y, t, dydt, part_a, part_b, relevant_disease);
        }
    }

    /**
     * @brief Checks whether LctState of a group satisfies all constraints. 
     *  Recursively, it checks that all groups satisfy the constraints.
     *
     * @tparam group The group for which the constraints should be checked.
     * @return Returns true if one or more constraints are not satisfied, false otherwise.
     */
    template <size_t Group = 0>
    bool check_constraints_impl() const
    {
        static_assert((Group < num_groups) && (Group >= 0), "The template parameter Group should be valid.");
        using LctStateGroup = type_at_index_t<Group, LctStates...>;

        if (LctStateGroup::template get_num_subcompartments<InfectionState::Susceptible>() != 1) {
            log_warning("Constraint check: The number of subcompartments for Susceptibles of group {} should be one!",
                        Group);
            return true;
        }
        if (LctStateGroup::template get_num_subcompartments<InfectionState::Recovered_1a>() != 1 ||
            LctStateGroup::template get_num_subcompartments<InfectionState::Recovered_1b>() != 1 ||
            LctStateGroup::template get_num_subcompartments<InfectionState::Recovered_2ab>() != 1) {
            log_warning("Constraint check: The number of subcompartments for Recovered of group {} should be one!",
                        Group);
            return true;
        }
        if (LctStateGroup::template get_num_subcompartments<InfectionState::Dead_a>() != 1 ||
            LctStateGroup::template get_num_subcompartments<InfectionState::Dead_b>() != 1) {
            log_warning("Constraint check: The number of subcompartments for Dead of group {} should be one!", Group);
            return true;
        }

        if constexpr (Group == num_groups - 1) {
            return false;
        }
        else {
            return check_constraints_impl<Group + 1>();
        }
    }
};

} // namespace lsecir2d
} // namespace mio
#endif // LCT_SECIR_2_DISEASE_MODEL_H
