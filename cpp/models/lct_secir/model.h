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

namespace mio
{
namespace lsecir
{

/**
 * @brief Class that defines an LCT-SECIR model.
 *
 * @tparam NumExposed The number of subcompartents used for the Exposed compartment.
 * @tparam NumInfectedNoSymptoms The number of subcompartents used for the InfectedNoSymptoms compartment. 
 * @tparam NumInfectedSymptoms The number of subcompartents used for the InfectedSymptoms compartment.
 * @tparam NumInfectedSevere The number of subcompartents used for the InfectedSevere compartment.
 * @tparam NumInfectedCritical The number of subcompartents used for the InfectedCritical compartment.
 */
template <typename FP = ScalarType, class... LctStates>
class Model
    : public CompartmentalModel<ScalarType, InfectionState, LctPopulations<ScalarType, LctStates...>, Parameters>
{ //TODO: Check constraints with check that num of subcompartments in S, R, D is one
public:
    using tupleLctStates = std::tuple<LctStates...>;
    using Base = CompartmentalModel<ScalarType, InfectionState, LctPopulations<ScalarType, LctStates...>, Parameters>;
    using typename Base::ParameterSet;
    using typename Base::Populations;
    static size_t constexpr m_elements = sizeof...(LctStates);
    assert(m_elements >= 1);

    /// @brief Default constructor.
    Model()
        : Base(Populations(), ParameterSet(m_elements))
    {
    }

    /// @brief TODO
    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    /**
     * @brief Evaluates the right-hand-side f of the LCT dydt = f(y, t).
     *
     * The LCT-SECIR model is defined through ordinary differential equations of the form dydt = f(y, t). 
     * y is a vector containing number of individuals for each (sub-) compartment.
     * This function evaluates the right-hand-side f of the ODE and can be used in an ODE solver.
     * @param[in] pop the current state of the population in the geographic unit we are considering
     * @param[in] y the current state of the model (or a subpopulation) as a flat array
     * @param[in] t the current time
     * @param[out] dydt a reference to the calculated output
     */
    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, ScalarType t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        // Vectors are sorted such that we first have all InfectionState%s for AgeGroup 0,
        // afterwards all for AgeGroup 1 and so on.
        dydt.setZero();
        get_derivatives_impl(pop, y, t, dydt);
    }
    template <size_t element1, size_t element2 = 0>
    ScalarType interact(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, ScalarType t,
                        Eigen::Ref<Eigen::VectorXd> dydt)
    {
        if constexpr ((element1 != element2) && (element2 < m_elements)) {
            size_t Si                       = this->populations.template get_first_index_element<element1>();
            ScalarType infectedNoSymptoms_j = 0;
            ScalarType infectedSymptoms_j   = 0;
            auto params                     = this->parameters;

            size_t elem2_first_index = this->populations.template get_first_index_element<element2>();

            // Calculate sum of all subcompartments for InfectedNoSymptoms of AgeGroup j.
            infectedNoSymptoms_j =
                pop.segment(elem2_first_index +
                                std::tuple_element_t<element2, tupleLctStates>::template get_first_index<
                                    InfectionState::InfectedNoSymptoms>(),
                            std::tuple_element_t<element2, tupleLctStates>::template get_num_subcompartments<
                                InfectionState::InfectedNoSymptoms>())
                    .sum();
            // Calculate sum of all subcompartments for InfectedSymptoms of AgeGroup j.
            infectedSymptoms_j =
                pop.segment(elem2_first_index +
                                std::tuple_element_t<element2, tupleLctStates>::template get_first_index<
                                    InfectionState::InfectedSymptoms>(),
                            std::tuple_element_t<element2, tupleLctStates>::template get_num_subcompartments<
                                InfectionState::InfectedSymptoms>())
                    .sum();
            // Size of the Subpopulation without dead people.
            double Nj = pop.segment(elem2_first_index, std::tuple_element_t<element2, tupleLctStates>::Count - 1).sum();
            ScalarType season_val =
                1 + params.template get<Seasonality>() *
                        sin(3.141592653589793 * ((params.template get<StartDay>() + t) / 182.5 + 0.5));
            dydt[Si] += -y[Si] / Nj * season_val * params.template get<TransmissionProbabilityOnContact>()[element1] *
                        params.template get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(
                            static_cast<Eigen::Index>(element1), static_cast<Eigen::Index>(element2)) *
                        (params.template get<RelativeTransmissionNoSymptoms>()[element2] * infectedNoSymptoms_j +
                         params.template get<RiskOfInfectionFromSymptomatic>()[element2] * infectedSymptoms_j);

            interact<element1, element2 + 1>(pop, y, t, dydt);
        }
    }
    template <size_t element = 0>
    void get_derivatives_impl(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, ScalarType t,
                              Eigen::Ref<Eigen::VectorXd> dydt) const
    {
        if constexpr (element > 0 && element < m_elements) {
            size_t first_index = this->populations.template get_first_index_element<element>();
            auto params        = this->parameters;
            ScalarType flow    = 0;

            // Indizes of first subcompartment of the InfectionState for the AgeGroup in an Vector.

            size_t Ei_first_index = this->populations.get_flat_index(
                {i, Index<LctState>(LctState::template get_first_index<InfectionState::Exposed>())});
            size_t INSi_first_index = this->populations.get_flat_index(
                {i, Index<LctState>(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>())});
            size_t ISyi_first_index = this->populations.get_flat_index(
                {i, Index<LctState>(LctState::template get_first_index<InfectionState::InfectedSymptoms>())});
            size_t ISevi_first_index = this->populations.get_flat_index(
                {i, Index<LctState>(LctState::template get_first_index<InfectionState::InfectedSevere>())});
            size_t ICri_first_index = this->populations.get_flat_index(
                {i, Index<LctState>(LctState::template get_first_index<InfectionState::InfectedCritical>())});
            size_t Ri = this->populations.get_flat_index(
                {i, Index<LctState>(LctState::template get_first_index<InfectionState::Recovered>())});
            size_t Di = this->populations.get_flat_index(
                {i, Index<LctState>(LctState::template get_first_index<InfectionState::Dead>())});

            // S'
            // AgeGroup i interacts with AgeGroup j.
            interact<element>();

            // E'
            dydt[Ei_first_index] = -dydt[Si];
            for (size_t subcomp = 0; subcomp < LctState::template get_num_subcompartments<InfectionState::Exposed>();
                 subcomp++) {
                // variable flow stores the value of the flow from one subcompartment to the next one.
                // Ei_first_index+subcomp is always the index of a (sub-)compartment of Exposed and Ei_first_index + subcomp+1 can also be the index of the first (sub-)compartment of InfectedNoSymptoms.
                flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::Exposed>() *
                       (1 / params.template get<TimeExposed>()[i]) * y[Ei_first_index + subcomp];
                // Subtract flow from dydt[1 + i] and add to dydt[2 + i].
                dydt[Ei_first_index + subcomp]     = dydt[Ei_first_index + subcomp] - flow;
                dydt[Ei_first_index + subcomp + 1] = flow;
            }

            // C'
            for (size_t subcomp = 0;
                 subcomp < LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>();
                 subcomp++) {
                flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>() *
                       (1 / params.template get<TimeInfectedNoSymptoms>()[i]) * y[INSi_first_index + subcomp];
                dydt[INSi_first_index + subcomp] -= flow;
                dydt[INSi_first_index + subcomp + 1] = flow;
            }

            // I'
            // Flow from last (sub-) compartment of C must be split between I_1 and R.
            dydt[Ri] = dydt[ISyi_first_index] * params.template get<RecoveredPerInfectedNoSymptoms>()[i];
            dydt[ISyi_first_index] =
                dydt[ISyi_first_index] * (1 - params.template get<RecoveredPerInfectedNoSymptoms>()[i]);

            for (size_t subcomp = 0;
                 subcomp < LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>(); subcomp++) {
                flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>() *
                       (1 / params.template get<TimeInfectedSymptoms>()[i]) * y[ISyi_first_index + subcomp];
                dydt[ISyi_first_index + subcomp] -= flow;
                dydt[ISyi_first_index + subcomp + 1] = flow;
            }

            // H'
            dydt[Ri] += dydt[ISevi_first_index] * (1 - params.template get<SeverePerInfectedSymptoms>()[i]);
            dydt[ISevi_first_index] = dydt[ISevi_first_index] * params.template get<SeverePerInfectedSymptoms>()[i];
            for (size_t subcomp = 0;
                 subcomp < LctState::template get_num_subcompartments<InfectionState::InfectedSevere>(); subcomp++) {
                flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSevere>() *
                       (1 / params.template get<TimeInfectedSevere>()[i]) * y[ISevi_first_index + subcomp];
                dydt[ISevi_first_index + subcomp] -= flow;
                dydt[ISevi_first_index + subcomp + 1] = flow;
            }

            // U'
            dydt[Ri] += dydt[ICri_first_index] * (1 - params.template get<CriticalPerSevere>()[i]);
            dydt[ICri_first_index] = dydt[ICri_first_index] * params.template get<CriticalPerSevere>()[i];
            for (size_t subcomp = 0;
                 subcomp < LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() - 1;
                 subcomp++) {
                flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() *
                       (1 / params.template get<TimeInfectedCritical>()[i]) * y[ICri_first_index + subcomp];
                dydt[ICri_first_index + subcomp] -= flow;
                dydt[ICri_first_index + subcomp + 1] = flow;
            }
            // Last flow from U has to be divided between R and D.
            // Must be calculated separately in order not to overwrite the already calculated values ​​for R.
            flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() *
                   (1 / params.template get<TimeInfectedCritical>()[i]) * y[Ri - 1];
            dydt[Ri - 1] -= flow;
            dydt[Ri] = dydt[Ri] + (1 - params.template get<DeathsPerCritical>()[i]) * flow;
            dydt[Di] = params.template get<DeathsPerCritical>()[i] * flow;

            get_derivatives_impl<element + 1>(pop, y, t, dydt)
        } // end if
    }

    /**
     * @brief Cumulates a simulation result with subcompartments to produce a result that divides the population only
     *   into the infection states defined in InfectionState.
     *
     * If the model is used for simulation, we will get a result in form of a TimeSeries with infection states divided 
     * in subcompartments.
     * The function calculates a TimeSeries without subcompartmens from another TimeSeries with subcompartments. 
     * This is done by summing up the numbers in the subcompartments.
     * @param[in] subcompartments_ts Result of a simulation with the model.
     * @return Result of the simulation divided in infection states without subcompartments. 
     *  Returns TimeSeries with values -1 if calculation is not possible.
     */
    // TimeSeries<ScalarType> calculate_compartments(const TimeSeries<ScalarType>& subcompartments_ts) const
    // {
    //     Eigen::Index count_InfStates = (Eigen::Index)InfectionState::Count;
    //     Eigen::Index num_compartments =
    //         count_InfStates * static_cast<Eigen::Index>((size_t)this->parameters.get_num_groups());
    //     TimeSeries<ScalarType> compartments_ts(num_compartments);
    //     if (!(this->populations.get_num_compartments() == (size_t)subcompartments_ts.get_num_elements())) {
    //         log_error("Result does not match infectionState of the Model.");
    //         Eigen::VectorXd wrong_size = Eigen::VectorXd::Constant(num_compartments, -1);
    //         compartments_ts.add_time_point(-1, wrong_size);
    //         return compartments_ts;
    //     }
    //     Eigen::VectorXd compartments(num_compartments);
    //     for (Eigen::Index timepoint = 0; timepoint < subcompartments_ts.get_num_time_points(); ++timepoint) {
    //         for (auto agegroup = AgeGroup(0); agegroup < this->parameters.get_num_groups(); agegroup++) {
    //             // Indizes of first subcompartment of the InfectionState for the AgeGroup in an Vector.
    //             size_t S_agegroup = this->populations.get_flat_index(
    //                 {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::Susceptible>())});
    //             size_t E_agegroup_first_index = this->populations.get_flat_index(
    //                 {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::Exposed>())});
    //             size_t INS_agegroup_first_index = this->populations.get_flat_index(
    //                 {agegroup,
    //                  Index<LctState>(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>())});
    //             size_t ISy_agegroup_first_index = this->populations.get_flat_index(
    //                 {agegroup,
    //                  Index<LctState>(LctState::template get_first_index<InfectionState::InfectedSymptoms>())});
    //             size_t ISev_agegroup_first_index = this->populations.get_flat_index(
    //                 {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::InfectedSevere>())});
    //             size_t ICr_agegroup_first_index = this->populations.get_flat_index(
    //                 {agegroup,
    //                  Index<LctState>(LctState::template get_first_index<InfectionState::InfectedCritical>())});
    //             size_t R_agegroup = this->populations.get_flat_index(
    //                 {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::Recovered>())});
    //             size_t D_agegroup = this->populations.get_flat_index(
    //                 {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::Dead>())});
    //             // Use segment of vector of the result with subcompartments of InfectionState with index j and sum up values of subcompartments.
    //             compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
    //                          (Eigen::Index)InfectionState::Susceptible] = subcompartments_ts[timepoint][S_agegroup];
    //             compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
    //                          (Eigen::Index)InfectionState::Exposed] =
    //                 subcompartments_ts[timepoint]
    //                     .segment(E_agegroup_first_index,
    //                              LctState::template get_num_subcompartments<InfectionState::Exposed>())
    //                     .sum();
    //             compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
    //                          (Eigen::Index)InfectionState::InfectedNoSymptoms] =
    //                 subcompartments_ts[timepoint]
    //                     .segment(INS_agegroup_first_index,
    //                              LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
    //                     .sum();
    //             compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
    //                          (Eigen::Index)InfectionState::InfectedSymptoms] =
    //                 subcompartments_ts[timepoint]
    //                     .segment(ISy_agegroup_first_index,
    //                              LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
    //                     .sum();
    //             compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
    //                          (Eigen::Index)InfectionState::InfectedSevere] =
    //                 subcompartments_ts[timepoint]
    //                     .segment(ISev_agegroup_first_index,
    //                              LctState::template get_num_subcompartments<InfectionState::InfectedSevere>())
    //                     .sum();
    //             compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
    //                          (Eigen::Index)InfectionState::InfectedCritical] =
    //                 subcompartments_ts[timepoint]
    //                     .segment(ICr_agegroup_first_index,
    //                              LctState::template get_num_subcompartments<InfectionState::InfectedCritical>())
    //                     .sum();
    //             compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
    //                          (Eigen::Index)InfectionState::Recovered] = subcompartments_ts[timepoint][R_agegroup];
    //             compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
    //                          (Eigen::Index)InfectionState::Dead]      = subcompartments_ts[timepoint][D_agegroup];
    //         } // end for agegroup

    //         compartments_ts.add_time_point(subcompartments_ts.get_time(timepoint), compartments);
    //     } // end for timepoints

    //     return compartments_ts;
    // }
};

} // namespace lsecir
} // namespace mio

#endif // LCTSECIR_MODEL_H
