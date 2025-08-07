/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Annika Jungkalus, Lena Ploetzke
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
#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/lct2d_populations.h"
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
template <class... LctStates>
class Model
    : public CompartmentalModel<ScalarType, InfectionState, LctPopulations<ScalarType, LctStates...>, Parameters>
{
public:
    using LctStatesGroups = TypeList<LctStates...>;
    using Base = CompartmentalModel<ScalarType, InfectionState, LctPopulations<ScalarType, LctStates...>, Parameters>;
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
    void get_derivatives(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop,
                         Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType t,
                         Eigen::Ref<Eigen::VectorX<ScalarType>> dydt) const override
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
    TimeSeries<ScalarType> calculate_compartments(const TimeSeries<ScalarType>& subcompartments_ts) const
    {
        Eigen::Index count_InfStates  = (Eigen::Index)InfectionState::Count;
        Eigen::Index num_compartments = count_InfStates * num_groups;
        TimeSeries<ScalarType> compartments_ts(num_compartments);
        if (!(this->populations.get_num_compartments() == (size_t)subcompartments_ts.get_num_elements())) {
            log_error("Result does not match InfectionState of the Model.");
            Eigen::VectorX<ScalarType> wrong_size = Eigen::VectorX<ScalarType>::Constant(num_compartments, -1);
            compartments_ts.add_time_point(-1, wrong_size);
            return compartments_ts;
        }
        Eigen::VectorX<ScalarType> compartments(num_compartments);
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
    void compress_vector(const Eigen::VectorX<ScalarType>& subcompartments,
                         Eigen::VectorX<ScalarType>& compartments) const
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
    void get_derivatives_impl(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop,
                              Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType t,
                              Eigen::Ref<Eigen::VectorX<ScalarType>> dydt) const
    {
        static_assert((Group < num_groups) && (Group >= 0), "The template parameter Group should be valid.");
        using LctStateGroup = type_at_index_t<Group, LctStates...>;

        size_t first_index_group = this->populations.template get_first_index_of_group<Group>();
        auto params              = this->parameters;
        ScalarType flow          = 0;

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
        size_t Ri_a  = first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered_a>();
        size_t Ri_b  = first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered_b>();
        size_t Ri_ab = first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered_ab>();
        size_t Di_a  = first_index_group + LctStateGroup::template get_first_index<InfectionState::Dead_a>();
        size_t Di_b  = first_index_group + LctStateGroup::template get_first_index<InfectionState::Dead_b>();

        // Calculate derivative of the Susceptible compartment.
        // outflow generated by disease a and disease b both
        double part_a = 0.;
        double part_b = 0.;
        interact<Group, 0>(pop, y, t, dydt, &part_a, &part_b, first_index_group, 2);
        // split flow
        double div_part_both =
            ((part_a + part_b) < Limits<ScalarType>::zero_tolerance()) ? 0.0 : 1.0 / (part_a + part_b);

        // Start with derivatives of First Infections (X_1a, X_1b)
        // Calculate derivative of the Exposed_1x compartments, split flow from Susceptible.
        // Exposed 1 a:
        dydt[Ei_1a_first_index] = -dydt[first_index_group] * part_a * div_part_both;
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_1a>(); subcomp++) {
            // Variable flow stores the value of the flow from one subcompartment to the next one.
            // Ei_1a_first_index + subcomp is always the index of a (sub-)compartment of Exposed and Ei_1a_first_index
            // + subcomp + 1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_1a>() *
                   (1 / params.template get<TimeExposed_a>()[Group]) * y[Ei_1a_first_index + subcomp];
            // Subtract flow from dydt[Ei_1a_first_index + subcomp] and add to next subcompartment.
            dydt[Ei_1a_first_index + subcomp] -= flow;
            dydt[Ei_1a_first_index + subcomp + 1] = flow;
        }
        // Exposed 1 b:
        dydt[Ei_1b_first_index] = -dydt[first_index_group] * part_b * div_part_both;
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_1b>(); subcomp++) {
            // Variable flow stores the value of the flow from one subcompartment to the next one.
            // Ei_1a_first_index + subcomp is always the index of a (sub-)compartment of Exposed and Ei_1b_first_index
            // + subcomp + 1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_1b>() *
                   (1 / params.template get<TimeExposed_b>()[Group]) * y[Ei_1b_first_index + subcomp];
            // Subtract flow from dydt[Ei_1b_first_index + subcomp] and add to next subcompartment.
            dydt[Ei_1b_first_index + subcomp] -= flow;
            dydt[Ei_1b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedNoSymptoms_1a compartment.
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1a>();
             subcomp++) {
            flow =
                (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1a>() *
                (1 / params.template get<TimeInfectedNoSymptoms_a>()[Group]) * y[INSi_1a_first_index + subcomp];
            dydt[INSi_1a_first_index + subcomp] -= flow;
            dydt[INSi_1a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedNoSymptoms_1b compartment.
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1b>();
             subcomp++) {
            flow =
                (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1b>() *
                (1 / params.template get<TimeInfectedNoSymptoms_b>()[Group]) * y[INSi_1b_first_index + subcomp];
            dydt[INSi_1b_first_index + subcomp] -= flow;
            dydt[INSi_1b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSymptoms_1a compartment.
        // Flow from last (sub-) compartment of C_1a must be split between
        // the first subcompartment of InfectedSymptoms_1a and Recovered_a.
        dydt[Ri_a] = dydt[ISyi_1a_first_index] * params.template get<RecoveredPerInfectedNoSymptoms_a>()[Group];
        dydt[ISyi_1a_first_index] =
            dydt[ISyi_1a_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms_a>()[Group]);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_1a>();
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_1a>() *
                   (1 / params.template get<TimeInfectedSymptoms_a>()[Group]) * y[ISyi_1a_first_index + subcomp];
            dydt[ISyi_1a_first_index + subcomp] -= flow;
            dydt[ISyi_1a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedSymptoms_1b compartment.
        // Flow from last (sub-) compartment of C_1b must be split between
        // the first subcompartment of InfectedSymptoms_1b and Recovered_b.
        dydt[Ri_b] = dydt[ISyi_1b_first_index] * params.template get<RecoveredPerInfectedNoSymptoms_b>()[Group];
        dydt[ISyi_1b_first_index] =
            dydt[ISyi_1b_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms_b>()[Group]);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_1b>();
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_1b>() *
                   (1 / params.template get<TimeInfectedSymptoms_b>()[Group]) * y[ISyi_1b_first_index + subcomp];
            dydt[ISyi_1b_first_index + subcomp] -= flow;
            dydt[ISyi_1b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSevere_1a compartment.
        // again split the flow from the last subcompartment of I_1a
        dydt[Ri_a] += dydt[ISevi_1a_first_index] * (1 - params.template get<SeverePerInfectedSymptoms_a>()[Group]);
        dydt[ISevi_1a_first_index] =
            dydt[ISevi_1a_first_index] * params.template get<SeverePerInfectedSymptoms_a>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_1a>();
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_1a>() *
                   (1 / params.template get<TimeInfectedSevere_a>()[Group]) * y[ISevi_1a_first_index + subcomp];
            dydt[ISevi_1a_first_index + subcomp] -= flow;
            dydt[ISevi_1a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedSevere_1b compartment.
        // again split the flow from the last subcompartment of I_1b
        dydt[Ri_b] += dydt[ISevi_1b_first_index] * (1 - params.template get<SeverePerInfectedSymptoms_b>()[Group]);
        dydt[ISevi_1b_first_index] =
            dydt[ISevi_1b_first_index] * params.template get<SeverePerInfectedSymptoms_b>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_1b>();
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_1b>() *
                   (1 / params.template get<TimeInfectedSevere_b>()[Group]) * y[ISevi_1b_first_index + subcomp];
            dydt[ISevi_1b_first_index + subcomp] -= flow;
            dydt[ISevi_1b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedCritical compartment.
        // again split flow from last subcompartment of H_1a between R_a and U_1a
        dydt[Ri_a] += dydt[ICri_1a_first_index] * (1 - params.template get<CriticalPerSevere_a>()[Group]);
        dydt[ICri_1a_first_index] = dydt[ICri_1a_first_index] * params.template get<CriticalPerSevere_a>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() - 1;
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() *
                   (1 / params.template get<TimeInfectedCritical_a>()[Group]) * y[ICri_1a_first_index + subcomp];
            dydt[ICri_1a_first_index + subcomp] -= flow;
            dydt[ICri_1a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedCritical compartment.
        // again split flow from last subcompartment of H_1b between R_b and U_1b
        dydt[Ri_b] += dydt[ICri_1b_first_index] * (1 - params.template get<CriticalPerSevere_b>()[Group]);
        dydt[ICri_1b_first_index] = dydt[ICri_1b_first_index] * params.template get<CriticalPerSevere_b>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() - 1;
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() *
                   (1 / params.template get<TimeInfectedCritical_b>()[Group]) * y[ICri_1b_first_index + subcomp];
            dydt[ICri_1b_first_index + subcomp] -= flow;
            dydt[ICri_1b_first_index + subcomp + 1] = flow;
        }

        // Last flow from InfectedCritical has to be divided between Recovered and Dead.
        // Must be calculated separately in order not to overwrite the already calculated values ​​for Recovered.
        // for 1a
        flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() *
               (1 / params.template get<TimeInfectedCritical_a>()[Group]) *
               y[ICri_1a_first_index +
                 LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() - 1];
        dydt[ICri_1a_first_index +
             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1a>() - 1] -= flow;
        dydt[Ri_a] = dydt[Ri_a] + (1 - params.template get<DeathsPerCritical_a>()[Group]) * flow;
        dydt[Di_a] = dydt[Di_a] + params.template get<DeathsPerCritical_a>()[Group] * flow;
        // for 1b
        flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() *
               (1 / params.template get<TimeInfectedCritical_b>()[Group]) *
               y[ICri_1b_first_index +
                 LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() - 1];
        dydt[ICri_1b_first_index +
             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_1b>() - 1] -= flow;
        dydt[Ri_b] = dydt[Ri_b] + (1 - params.template get<DeathsPerCritical_b>()[Group]) * flow;
        dydt[Di_b] = dydt[Di_b] + params.template get<DeathsPerCritical_b>()[Group] * flow;

        // Second Infection (X_2a, X_2b)
        // outflow from Recovered, people getting infected for the second time
        double temp_Ra = dydt[Ri_a];
        interact<0, 0>(pop, y, t, dydt, &part_a, &part_b, Ri_a,
                       1); // outflow from R_a is only affected by b, no age groups
        double temp_Rb = dydt[Ri_b];
        interact<0, 0>(pop, y, t, dydt, &part_a, &part_b, Ri_b,
                       0); // outflow from R_b is only affected by a, no age groups

        // Calculate derivative of the Exposed_2i compartments
        // Exposed 2 a
        dydt[Ei_2a_first_index] = -(dydt[Ri_b] - temp_Rb);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_2a>(); subcomp++) {
            // Variable flow stores the value of the flow from one subcompartment to the next one.
            // Ei_2a_first_index + subcomp is always the index of a (sub-)compartment of Exposed and Ei_2a_first_index
            // + subcomp + 1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_2a>() *
                   (1 / params.template get<TimeExposed_a>()[Group]) * y[Ei_2a_first_index + subcomp];
            // Subtract flow from dydt[Ei_2a_first_index + subcomp] and add to next subcompartment.
            dydt[Ei_2a_first_index + subcomp] -= flow;
            dydt[Ei_2a_first_index + subcomp + 1] = flow;
        }
        // Exposed 2 b
        dydt[Ei_2b_first_index] = -(dydt[Ri_a] - temp_Ra);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_2b>(); subcomp++) {
            // Variable flow stores the value of the flow from one subcompartment to the next one.
            // Ei_2b_first_index + subcomp is always the index of a (sub-)compartment of Exposed and Ei_2b_first_index
            // + subcomp + 1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed_2b>() *
                   (1 / params.template get<TimeExposed_b>()[Group]) * y[Ei_2b_first_index + subcomp];
            // Subtract flow from dydt[Ei_2b_first_index + subcomp] and add to next subcompartment.
            dydt[Ei_2b_first_index + subcomp] -= flow;
            dydt[Ei_2b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedNoSymptoms_2a compartment.
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2a>();
             subcomp++) {
            flow =
                (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2a>() *
                (1 / params.template get<TimeInfectedNoSymptoms_a>()[Group]) * y[INSi_2a_first_index + subcomp];
            dydt[INSi_2a_first_index + subcomp] -= flow;
            dydt[INSi_2a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedNoSymptoms_2b compartment.
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2b>();
             subcomp++) {
            flow =
                (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2b>() *
                (1 / params.template get<TimeInfectedNoSymptoms_b>()[Group]) * y[INSi_2b_first_index + subcomp];
            dydt[INSi_2b_first_index + subcomp] -= flow;
            dydt[INSi_2b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSymptoms_2a compartment.
        // Flow from last (sub-) compartment of C_2a must be split between
        // the first subcompartment of InfectedSymptoms_2a and Recovered_ab.
        dydt[Ri_ab] += dydt[ISyi_2a_first_index] * params.template get<RecoveredPerInfectedNoSymptoms_a>()[Group];
        dydt[ISyi_2a_first_index] =
            dydt[ISyi_2a_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms_a>()[Group]);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_2a>();
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_2a>() *
                   (1 / params.template get<TimeInfectedSymptoms_a>()[Group]) * y[ISyi_2a_first_index + subcomp];
            dydt[ISyi_2a_first_index + subcomp] -= flow;
            dydt[ISyi_2a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedSymptoms_2b compartment.
        // Flow from last (sub-) compartment of C_2b must be split between
        // the first subcompartment of InfectedSymptoms_2b and Recovered_ab.
        dydt[Ri_ab] += dydt[ISyi_2b_first_index] * params.template get<RecoveredPerInfectedNoSymptoms_b>()[Group];
        dydt[ISyi_2b_first_index] =
            dydt[ISyi_2b_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms_b>()[Group]);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_2b>();
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms_2b>() *
                   (1 / params.template get<TimeInfectedSymptoms_b>()[Group]) * y[ISyi_2b_first_index + subcomp];
            dydt[ISyi_2b_first_index + subcomp] -= flow;
            dydt[ISyi_2b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSevere_2a compartment.
        // again split the flow from the last subcompartment of I_2a
        dydt[Ri_ab] += dydt[ISevi_2a_first_index] * (1 - params.template get<SeverePerInfectedSymptoms_a>()[Group]);
        dydt[ISevi_2a_first_index] =
            dydt[ISevi_2a_first_index] * params.template get<SeverePerInfectedSymptoms_a>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_2a>();
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_2a>() *
                   (1 / params.template get<TimeInfectedSevere_a>()[Group]) * y[ISevi_2a_first_index + subcomp];
            dydt[ISevi_2a_first_index + subcomp] -= flow;
            dydt[ISevi_2a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedSevere compartment.
        // again split the flow from the last subcompartment of I_2b
        dydt[Ri_ab] += dydt[ISevi_2b_first_index] * (1 - params.template get<SeverePerInfectedSymptoms_b>()[Group]);
        dydt[ISevi_2b_first_index] =
            dydt[ISevi_2b_first_index] * params.template get<SeverePerInfectedSymptoms_b>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_2b>();
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere_2b>() *
                   (1 / params.template get<TimeInfectedSevere_b>()[Group]) * y[ISevi_2b_first_index + subcomp];
            dydt[ISevi_2b_first_index + subcomp] -= flow;
            dydt[ISevi_2b_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedCritical compartment.
        // again split flow from last subcompartment of H_2a between R_ab and U_2a
        dydt[Ri_ab] += dydt[ICri_2a_first_index] * (1 - params.template get<CriticalPerSevere_a>()[Group]);
        dydt[ICri_2a_first_index] = dydt[ICri_2a_first_index] * params.template get<CriticalPerSevere_a>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() - 1;
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() *
                   (1 / params.template get<TimeInfectedCritical_a>()[Group]) * y[ICri_2a_first_index + subcomp];
            dydt[ICri_2a_first_index + subcomp] -= flow;
            dydt[ICri_2a_first_index + subcomp + 1] = flow;
        }
        // Calculate derivative of the InfectedCritical compartment.
        // again split flow from last subcompartment of H_1b between R_b and U_1b
        dydt[Ri_ab] += dydt[ICri_2b_first_index] * (1 - params.template get<CriticalPerSevere_b>()[Group]);
        dydt[ICri_2b_first_index] = dydt[ICri_2b_first_index] * params.template get<CriticalPerSevere_b>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() - 1;
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() *
                   (1 / params.template get<TimeInfectedCritical_b>()[Group]) * y[ICri_2b_first_index + subcomp];
            dydt[ICri_2b_first_index + subcomp] -= flow;
            dydt[ICri_2b_first_index + subcomp + 1] = flow;
        }

        // Last flow from InfectedCritical has to be divided between Recovered and Dead.
        // Must be calculated separately in order not to overwrite the already calculated values ​​for Recovered.
        // for 2a
        flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() *
               (1 / params.template get<TimeInfectedCritical_a>()[Group]) *
               y[ICri_2a_first_index +
                 LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() - 1];
        dydt[ICri_2a_first_index +
             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2a>() - 1] -= flow;
        dydt[Ri_ab] = dydt[Ri_ab] + (1 - params.template get<DeathsPerCritical_a>()[Group]) * flow;
        dydt[Di_a]  = dydt[Di_a] + params.template get<DeathsPerCritical_a>()[Group] * flow;
        // for 2b
        flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() *
               (1 / params.template get<TimeInfectedCritical_b>()[Group]) *
               y[ICri_2b_first_index +
                 LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() - 1];
        dydt[ICri_2b_first_index +
             LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical_2b>() - 1] -= flow;
        dydt[Ri_ab] = dydt[Ri_ab] + (1 - params.template get<DeathsPerCritical_b>()[Group]) * flow;
        dydt[Di_b]  = dydt[Di_b] + params.template get<DeathsPerCritical_b>()[Group] * flow;

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
     * @param[in] compartment_index Index of the compartment for that the outflow will be calculated.
     * @param[in] which_disease Index for which infected people cause the outflow (0 = a, 1 = b, 2 = a and b).
     */
    template <size_t Group1, size_t Group2 = 0>
    void interact(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                  ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> dydt, double* part_a, double* part_b,
                  size_t compartment_index, int which_disease) const
    {
        static_assert((Group1 < num_groups) && (Group1 >= 0) && (Group2 < num_groups) && (Group2 >= 0),
                      "The template parameters Group1 & Group2 should be valid.");
        using LctStateGroup2              = type_at_index_t<Group2, LctStates...>;
        ScalarType infectedNoSymptoms_2_a = 0;
        ScalarType infectedSymptoms_2_a   = 0;
        ScalarType infectedNoSymptoms_2_b = 0;
        ScalarType infectedSymptoms_2_b   = 0;
        auto params                       = this->parameters;

        size_t first_index_group2 = this->populations.template get_first_index_of_group<Group2>();

        // Calculate sum of all subcompartments for InfectedNoSymptoms for disease a of Group2.
        infectedNoSymptoms_2_a =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedNoSymptoms_1a>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1a>())
                .sum() +
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedNoSymptoms_2a>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2a>())
                .sum();
        // Calculate sum of all subcompartments for InfectedSymptoms for disease a of Group2.
        infectedSymptoms_2_a =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedSymptoms_1a>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedSymptoms_1a>())
                .sum() +
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedSymptoms_2a>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedSymptoms_2a>())
                .sum();
        // Calculate sum of all subcompartments for InfectedNoSymptoms for disease b of Group2.
        infectedNoSymptoms_2_b =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedNoSymptoms_1b>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_1b>())
                .sum() +
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedNoSymptoms_2b>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedNoSymptoms_2b>())
                .sum();
        // Calculate sum of all subcompartments for InfectedSymptoms for disease b of Group2.
        infectedSymptoms_2_b =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedSymptoms_1b>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedSymptoms_1b>())
                .sum() +
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedSymptoms_2b>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedSymptoms_2b>())
                .sum();
        // Size of the Subpopulation Group2 without dead people.
        double N_2 = pop.segment(first_index_group2, LctStateGroup2::Count).sum(); // sum over all compartments
        N_2        = N_2 - pop.segment(LctStateGroup2::template get_first_index<InfectionState::Dead_a>(), 1).sum() -
              pop.segment(LctStateGroup2::template get_first_index<InfectionState::Dead_b>(), 1)
                  .sum(); // all people minus dead people
        const double divN_2   = (N_2 < Limits<ScalarType>::zero_tolerance()) ? 0.0 : 1.0 / N_2;
        ScalarType season_val = 1 + params.template get<Seasonality>() *
                                        sin(3.141592653589793 * ((params.template get<StartDay>() + t) / 182.5 + 0.5));

        if (which_disease == 0) { // disease a
            dydt[compartment_index] +=
                -y[compartment_index] * divN_2 * season_val *
                params.template get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(
                    static_cast<Eigen::Index>(Group1), static_cast<Eigen::Index>(Group2)) *
                params.template get<TransmissionProbabilityOnContact_a>()[Group1] *
                (params.template get<RelativeTransmissionNoSymptoms_a>()[Group2] * infectedNoSymptoms_2_a +
                 params.template get<RiskOfInfectionFromSymptomatic_a>()[Group2] * infectedSymptoms_2_a);
        }
        else if (which_disease == 1) { // disease b
            dydt[compartment_index] +=
                -y[compartment_index] * divN_2 * season_val *
                params.template get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(
                    static_cast<Eigen::Index>(Group1), static_cast<Eigen::Index>(Group2)) *
                params.template get<TransmissionProbabilityOnContact_b>()[Group1] *
                (params.template get<RelativeTransmissionNoSymptoms_b>()[Group2] * infectedNoSymptoms_2_b +
                 params.template get<RiskOfInfectionFromSymptomatic_b>()[Group2] * infectedSymptoms_2_b);
        }
        else if (which_disease == 2) { // both diseases
            dydt[compartment_index] +=
                -y[compartment_index] * divN_2 * season_val *
                params.template get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(
                    static_cast<Eigen::Index>(Group1), static_cast<Eigen::Index>(Group2)) *
                (params.template get<TransmissionProbabilityOnContact_a>()[Group1] *
                     (params.template get<RelativeTransmissionNoSymptoms_a>()[Group2] * infectedNoSymptoms_2_a +
                      params.template get<RiskOfInfectionFromSymptomatic_a>()[Group2] * infectedSymptoms_2_a) +
                 params.template get<TransmissionProbabilityOnContact_b>()[Group1] *
                     (params.template get<RelativeTransmissionNoSymptoms_b>()[Group2] * infectedNoSymptoms_2_b +
                      params.template get<RiskOfInfectionFromSymptomatic_b>()[Group2] * infectedSymptoms_2_b));
        }
        // To split the outflow from S between E_1a and E_1b:
        *part_a += params.template get<TransmissionProbabilityOnContact_a>()[Group1] *
                   (params.template get<RelativeTransmissionNoSymptoms_a>()[Group2] * infectedNoSymptoms_2_a +
                    params.template get<RiskOfInfectionFromSymptomatic_a>()[Group2] * infectedSymptoms_2_a);
        *part_b += params.template get<TransmissionProbabilityOnContact_b>()[Group1] *
                   (params.template get<RelativeTransmissionNoSymptoms_b>()[Group2] * infectedNoSymptoms_2_b +
                    params.template get<RiskOfInfectionFromSymptomatic_b>()[Group2] * infectedSymptoms_2_b);

        if constexpr (Group2 + 1 < num_groups) {
            interact<Group1, Group2 + 1>(pop, y, t, dydt, part_a, part_b, compartment_index, which_disease);
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
        if (LctStateGroup::template get_num_subcompartments<InfectionState::Recovered_a>() != 1 ||
            LctStateGroup::template get_num_subcompartments<InfectionState::Recovered_b>() != 1 ||
            LctStateGroup::template get_num_subcompartments<InfectionState::Recovered_ab>() != 1) {
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
