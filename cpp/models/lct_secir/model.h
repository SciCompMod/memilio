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
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/age_group.h"
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
template <size_t NumExposed, size_t NumInfectedNoSymptoms, size_t NumInfectedSymptoms, size_t NumInfectedSevere,
          size_t NumInfectedCritical>
class Model
    : public CompartmentalModel<
          ScalarType,
          LctInfectionState<InfectionState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                            NumInfectedSevere, NumInfectedCritical, 1, 1>,
          mio::Populations<ScalarType, AgeGroup,
                           LctInfectionState<InfectionState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                                             NumInfectedSevere, NumInfectedCritical, 1, 1>>,
          Parameters>
{
public:
    using LctState = LctInfectionState<InfectionState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                                       NumInfectedSevere, NumInfectedCritical, 1, 1>;
    using Base = CompartmentalModel<ScalarType, LctState, mio::Populations<ScalarType, AgeGroup, LctState>, Parameters>;
    using typename Base::ParameterSet;
    using typename Base::Populations;

    /// @brief Default constructor.
    //TODO: why do we need the 0.0 here?
    Model()
        : Base(Populations({AgeGroup(1), Index<LctState>(LctState::Count)}, 0.), ParameterSet(AgeGroup(1)))
    {
    }

    /// @brief TODO
    Model(int num_agegroups)
        : Base(Populations({AgeGroup(num_agegroups), Index<LctState>(LctState::Count)}, 0.),
               ParameterSet(AgeGroup(num_agegroups)))
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

        auto params                     = this->parameters;
        ScalarType infectedNoSymptoms_j = 0;
        ScalarType infectedSymptoms_j   = 0;
        ScalarType flow                 = 0;

        for (auto i = AgeGroup(0); i < params.get_num_groups(); i++) {
            // Indizes of first subcompartment of the InfectionState for the AgeGroup in an Vector.
            size_t Si = this->populations.get_flat_index(
                {i, Index<LctState>(LctState::template get_first_index<InfectionState::Susceptible>())});
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
            for (auto j = AgeGroup(0); j < params.get_num_groups(); j++) {
                size_t INSj_first_index = this->populations.get_flat_index(
                    {j, Index<LctState>(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>())});
                size_t ISyj_first_index = this->populations.get_flat_index(
                    {j, Index<LctState>(LctState::template get_first_index<InfectionState::InfectedSymptoms>())});
                size_t agegroupj_first_index = this->populations.get_flat_index(
                    {j, Index<LctState>(LctState::template get_first_index<InfectionState::Susceptible>())});
                // Calculate sum of all subcompartments for InfectedNoSymptoms of AgeGroup j.
                infectedNoSymptoms_j =
                    pop.segment(INSj_first_index,
                                LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
                        .sum();
                // Calculate sum of all subcompartments for InfectedSymptoms of AgeGroup j.
                infectedSymptoms_j =
                    pop.segment(ISyj_first_index,
                                LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
                        .sum();
                // Size of the Subpopulation without dead people.
                double Nj = pop.segment(agegroupj_first_index, LctState::Count - 1).sum();
                ScalarType season_val =
                    1 + params.template get<Seasonality>() *
                            sin(3.141592653589793 * ((params.template get<StartDay>() + t) / 182.5 + 0.5));
                dydt[Si] += -y[Si] / Nj * season_val * params.template get<TransmissionProbabilityOnContact>()[i] *
                            params.template get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(
                                static_cast<Eigen::Index>((size_t)i), static_cast<Eigen::Index>((size_t)j)) *
                            (params.template get<RelativeTransmissionNoSymptoms>()[j] * infectedNoSymptoms_j +
                             params.template get<RiskOfInfectionFromSymptomatic>()[j] * infectedSymptoms_j);

                infectedNoSymptoms_j = 0;
                infectedSymptoms_j   = 0;
            }

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
        } // end for
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
    TimeSeries<ScalarType> calculate_compartments(const TimeSeries<ScalarType>& subcompartments_ts) const
    {
        Eigen::Index count_InfStates = (Eigen::Index)InfectionState::Count;
        Eigen::Index num_compartments =
            count_InfStates * static_cast<Eigen::Index>((size_t)this->parameters.get_num_groups());
        TimeSeries<ScalarType> compartments_ts(num_compartments);
        if (!(this->populations.get_num_compartments() == (size_t)subcompartments_ts.get_num_elements())) {
            log_error("Result does not match infectionState of the Model.");
            Eigen::VectorXd wrong_size = Eigen::VectorXd::Constant(num_compartments, -1);
            compartments_ts.add_time_point(-1, wrong_size);
            return compartments_ts;
        }
        Eigen::VectorXd compartments(num_compartments);
        for (Eigen::Index timepoint = 0; timepoint < subcompartments_ts.get_num_time_points(); ++timepoint) {
            for (auto agegroup = AgeGroup(0); agegroup < this->parameters.get_num_groups(); agegroup++) {
                // Indizes of first subcompartment of the InfectionState for the AgeGroup in an Vector.
                size_t S_agegroup = this->populations.get_flat_index(
                    {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::Susceptible>())});
                size_t E_agegroup_first_index = this->populations.get_flat_index(
                    {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::Exposed>())});
                size_t INS_agegroup_first_index = this->populations.get_flat_index(
                    {agegroup,
                     Index<LctState>(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>())});
                size_t ISy_agegroup_first_index = this->populations.get_flat_index(
                    {agegroup,
                     Index<LctState>(LctState::template get_first_index<InfectionState::InfectedSymptoms>())});
                size_t ISev_agegroup_first_index = this->populations.get_flat_index(
                    {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::InfectedSevere>())});
                size_t ICr_agegroup_first_index = this->populations.get_flat_index(
                    {agegroup,
                     Index<LctState>(LctState::template get_first_index<InfectionState::InfectedCritical>())});
                size_t R_agegroup = this->populations.get_flat_index(
                    {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::Recovered>())});
                size_t D_agegroup = this->populations.get_flat_index(
                    {agegroup, Index<LctState>(LctState::template get_first_index<InfectionState::Dead>())});
                // Use segment of vector of the result with subcompartments of InfectionState with index j and sum up values of subcompartments.
                compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
                             (Eigen::Index)InfectionState::Susceptible] = subcompartments_ts[timepoint][S_agegroup];
                compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
                             (Eigen::Index)InfectionState::Exposed] =
                    subcompartments_ts[timepoint]
                        .segment(E_agegroup_first_index,
                                 LctState::template get_num_subcompartments<InfectionState::Exposed>())
                        .sum();
                compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
                             (Eigen::Index)InfectionState::InfectedNoSymptoms] =
                    subcompartments_ts[timepoint]
                        .segment(INS_agegroup_first_index,
                                 LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
                        .sum();
                compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
                             (Eigen::Index)InfectionState::InfectedSymptoms] =
                    subcompartments_ts[timepoint]
                        .segment(ISy_agegroup_first_index,
                                 LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
                        .sum();
                compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
                             (Eigen::Index)InfectionState::InfectedSevere] =
                    subcompartments_ts[timepoint]
                        .segment(ISev_agegroup_first_index,
                                 LctState::template get_num_subcompartments<InfectionState::InfectedSevere>())
                        .sum();
                compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
                             (Eigen::Index)InfectionState::InfectedCritical] =
                    subcompartments_ts[timepoint]
                        .segment(ICr_agegroup_first_index,
                                 LctState::template get_num_subcompartments<InfectionState::InfectedCritical>())
                        .sum();
                compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
                             (Eigen::Index)InfectionState::Recovered] = subcompartments_ts[timepoint][R_agegroup];
                compartments[count_InfStates * static_cast<Eigen::Index>((size_t)agegroup) +
                             (Eigen::Index)InfectionState::Dead]      = subcompartments_ts[timepoint][D_agegroup];
            } // end for agegroup

            compartments_ts.add_time_point(subcompartments_ts.get_time(timepoint), compartments);
        } // end for timepoints

        return compartments_ts;
    }
};

} // namespace lsecir
} // namespace mio

#endif // LCTSECIR_MODEL_H
