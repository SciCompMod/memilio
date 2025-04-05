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

#ifndef LCT_SECIR_MODEL_H
#define LCT_SECIR_MODEL_H

#include "lct_secir/parameters.h"
#include "lct_secir/infection_state.h"
#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/lct_populations.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/type_list.h"
#include "memilio/utils/metaprogramming.h"
namespace mio
{
namespace lsecir
{
/**
 * @brief Class that defines an LCT-SECIR model.
 *
 * @tparam LctStates The LCT model can work with any number of LctStates, where each LctState corresponds to a group,
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
     * The LCT-SECIR model is defined through ordinary differential equations of the form dydt = f(y, t). 
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
        size_t Ei_first_index = first_index_group + LctStateGroup::template get_first_index<InfectionState::Exposed>();
        size_t INSi_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedNoSymptoms>();
        size_t ISyi_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSymptoms>();
        size_t ISevi_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedSevere>();
        size_t ICri_first_index =
            first_index_group + LctStateGroup::template get_first_index<InfectionState::InfectedCritical>();
        size_t Ri = first_index_group + LctStateGroup::template get_first_index<InfectionState::Recovered>();
        size_t Di = first_index_group + LctStateGroup::template get_first_index<InfectionState::Dead>();

        // Calculate derivative of the Susceptible compartment.
        interact<Group, 0>(pop, y, t, dydt);

        // Calculate derivative of the Exposed compartment.
        dydt[Ei_first_index] = -dydt[first_index_group];
        for (size_t subcomp = 0; subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::Exposed>();
             subcomp++) {
            // Variable flow stores the value of the flow from one subcompartment to the next one.
            // Ei_first_index + subcomp is always the index of a (sub-)compartment of Exposed and Ei_first_index
            // + subcomp + 1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::Exposed>() *
                   (1 / params.template get<TimeExposed>()[Group]) * y[Ei_first_index + subcomp];
            // Subtract flow from dydt[Ei_first_index + subcomp] and add to next subcompartment.
            dydt[Ei_first_index + subcomp] -= flow;
            dydt[Ei_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedNoSymptoms compartment.
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>();
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>() *
                   (1 / params.template get<TimeInfectedNoSymptoms>()[Group]) * y[INSi_first_index + subcomp];
            dydt[INSi_first_index + subcomp] -= flow;
            dydt[INSi_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSymptoms compartment.
        // Flow from last (sub-) compartment of C must be split between
        // the first subcompartment of InfectedSymptoms and Recovered.
        dydt[Ri] = dydt[ISyi_first_index] * params.template get<RecoveredPerInfectedNoSymptoms>()[Group];
        dydt[ISyi_first_index] =
            dydt[ISyi_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms>()[Group]);
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms>(); subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSymptoms>() *
                   (1 / params.template get<TimeInfectedSymptoms>()[Group]) * y[ISyi_first_index + subcomp];
            dydt[ISyi_first_index + subcomp] -= flow;
            dydt[ISyi_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedSevere compartment.
        dydt[Ri] += dydt[ISevi_first_index] * (1 - params.template get<SeverePerInfectedSymptoms>()[Group]);
        dydt[ISevi_first_index] = dydt[ISevi_first_index] * params.template get<SeverePerInfectedSymptoms>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere>(); subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedSevere>() *
                   (1 / params.template get<TimeInfectedSevere>()[Group]) * y[ISevi_first_index + subcomp];
            dydt[ISevi_first_index + subcomp] -= flow;
            dydt[ISevi_first_index + subcomp + 1] = flow;
        }

        // Calculate derivative of the InfectedCritical compartment.
        dydt[Ri] += dydt[ICri_first_index] * (1 - params.template get<CriticalPerSevere>()[Group]);
        dydt[ICri_first_index] = dydt[ICri_first_index] * params.template get<CriticalPerSevere>()[Group];
        for (size_t subcomp = 0;
             subcomp < LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>() - 1;
             subcomp++) {
            flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>() *
                   (1 / params.template get<TimeInfectedCritical>()[Group]) * y[ICri_first_index + subcomp];
            dydt[ICri_first_index + subcomp] -= flow;
            dydt[ICri_first_index + subcomp + 1] = flow;
        }
        // Last flow from InfectedCritical has to be divided between Recovered and Dead.
        // Must be calculated separately in order not to overwrite the already calculated values ​​for Recovered.
        flow = (ScalarType)LctStateGroup::template get_num_subcompartments<InfectionState::InfectedCritical>() *
               (1 / params.template get<TimeInfectedCritical>()[Group]) * y[Ri - 1];
        dydt[Ri - 1] -= flow;
        dydt[Ri] = dydt[Ri] + (1 - params.template get<DeathsPerCritical>()[Group]) * flow;
        dydt[Di] = params.template get<DeathsPerCritical>()[Group] * flow;

        // Function call for next group if applicable.
        if constexpr (Group + 1 < num_groups) {
            get_derivatives_impl<Group + 1>(pop, y, t, dydt);
        }
    }

    /**
     * @brief Calculates the derivative of the Susceptible compartment for Group1.
     *
     * This is done recursively by calculating the interaction terms with each group.
     * @tparam Group1 The group for which the derivative of the Susceptible compartment should be calculated. 
     * @tparam Group2 The group that Group1 interacts with. 
     * @param[in] pop The current state of the population in the geographic unit we are considering.
     * @param[in] y The current state of the model (or a subpopulation) as a flat array.
     * @param[in] t The current time.
     * @param[out] dydt A reference to the calculated output.
     */
    template <size_t Group1, size_t Group2 = 0>
    void interact(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop, Eigen::Ref<const Eigen::VectorX<ScalarType>> y,
                  ScalarType t, Eigen::Ref<Eigen::VectorX<ScalarType>> dydt) const
    {
        static_assert((Group1 < num_groups) && (Group1 >= 0) && (Group2 < num_groups) && (Group2 >= 0),
                      "The template parameters Group1 & Group2 should be valid.");
        using LctStateGroup2            = type_at_index_t<Group2, LctStates...>;
        size_t Si_1                     = this->populations.template get_first_index_of_group<Group1>();
        ScalarType infectedNoSymptoms_2 = 0;
        ScalarType infectedSymptoms_2   = 0;
        auto params                     = this->parameters;

        size_t first_index_group2 = this->populations.template get_first_index_of_group<Group2>();

        // Calculate sum of all subcompartments for InfectedNoSymptoms of Group2.
        infectedNoSymptoms_2 =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
                .sum();
        // Calculate sum of all subcompartments for InfectedSymptoms of Group2.
        infectedSymptoms_2 =
            pop.segment(first_index_group2 +
                            LctStateGroup2::template get_first_index<InfectionState::InfectedSymptoms>(),
                        LctStateGroup2::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
                .sum();
        // Size of the Subpopulation Group2 without dead people.
        double N_2            = pop.segment(first_index_group2, LctStateGroup2::Count - 1).sum();
        const double divN_2   = (N_2 < Limits<ScalarType>::zero_tolerance()) ? 0.0 : 1.0 / N_2;
        ScalarType season_val = 1 + params.template get<Seasonality>() *
                                        sin(3.141592653589793 * ((params.template get<StartDay>() + t) / 182.5 + 0.5));
        dydt[Si_1] += -y[Si_1] * divN_2 * season_val * params.template get<TransmissionProbabilityOnContact>()[Group1] *
                      params.template get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(
                          static_cast<Eigen::Index>(Group1), static_cast<Eigen::Index>(Group2)) *
                      (params.template get<RelativeTransmissionNoSymptoms>()[Group2] * infectedNoSymptoms_2 +
                       params.template get<RiskOfInfectionFromSymptomatic>()[Group2] * infectedSymptoms_2);
        // Function call for next interacting group if applicable.
        if constexpr (Group2 + 1 < num_groups) {
            interact<Group1, Group2 + 1>(pop, y, t, dydt);
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
        if (LctStateGroup::template get_num_subcompartments<InfectionState::Recovered>() != 1) {
            log_warning("Constraint check: The number of subcompartments for Recovered of group {} should be one!",
                        Group);
            return true;
        }
        if (LctStateGroup::template get_num_subcompartments<InfectionState::Dead>() != 1) {
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

} // namespace lsecir
} // namespace mio

#endif // LCTSECIR_MODEL_H
