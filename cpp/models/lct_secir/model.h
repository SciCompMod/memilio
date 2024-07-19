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
          mio::Populations<ScalarType,
                           LctInfectionState<InfectionState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                                             NumInfectedSevere, NumInfectedCritical, 1, 1>>,
          Parameters>
{
public:
    using LctState = LctInfectionState<InfectionState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                                       NumInfectedSevere, NumInfectedCritical, 1, 1>;
    using Base     = CompartmentalModel<ScalarType, LctState, mio::Populations<ScalarType, LctState>, Parameters>;
    using typename Base::ParameterSet;
    using typename Base::Populations;

    /// @brief Default constructor.
    Model()
        : Base(Populations({Index<LctState>(LctState::Count)}, 0.), ParameterSet())
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
        dydt.setZero();

        auto params           = this->parameters;
        auto total_population = pop.sum() - pop[LctState::template get_first_index<InfectionState::Dead>()];

        ScalarType infectedNoSymptoms = 0;
        ScalarType infectedSymptoms   = 0;
        ScalarType flow               = 0;

        // Calculate sum of all subcompartments for InfectedNoSymptoms.
        infectedNoSymptoms =
            pop.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                        LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
                .sum();
        // Calculate sum of all subcompartments for InfectedSymptoms.
        infectedSymptoms = pop.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                                       LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
                               .sum();

        // S'
        ScalarType season_val = 1 + params.template get<Seasonality>() *
                                        sin(3.141592653589793 * ((params.template get<StartDay>() + t) / 182.5 + 0.5));
        dydt[0] = -y[0] / total_population * season_val * params.template get<TransmissionProbabilityOnContact>() *
                  params.template get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(0, 0) *
                  (params.template get<RelativeTransmissionNoSymptoms>() * infectedNoSymptoms +
                   params.template get<RiskOfInfectionFromSymptomatic>() * infectedSymptoms);

        // E'
        dydt[1] = -dydt[0];
        for (size_t i = 0; i < LctState::template get_num_subcompartments<InfectionState::Exposed>(); i++) {
            // Dummy stores the value of the flow from dydt[1 + i] to dydt[2 + i].
            // 1+i is always the index of a (sub-)compartment of E and 2+i can also be the index of the first (sub-)compartment of C.
            flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::Exposed>() *
                   (1 / params.template get<TimeExposed>()) * y[1 + i];
            // Subtract flow from dydt[1 + i] and add to dydt[2 + i].
            dydt[1 + i] = dydt[1 + i] - flow;
            dydt[2 + i] = flow;
        }

        // C'
        for (size_t i = 0; i < LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>(); i++) {
            flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>() *
                   (1 / params.template get<TimeInfectedNoSymptoms>()) *
                   y[LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() + i];
            dydt[LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() + i] =
                dydt[LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() + i] - flow;
            dydt[LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() + i + 1] = flow;
        }

        // I'
        // Flow from last (sub-) compartment of C must be split between I_1 and R.
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] =
            dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>()] *
            params.template get<RecoveredPerInfectedNoSymptoms>();
        dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>()] =
            dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>()] *
            (1 - params.template get<RecoveredPerInfectedNoSymptoms>());

        for (size_t i = 0; i < LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>(); i++) {
            flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>() *
                   (1 / params.template get<TimeInfectedSymptoms>()) *
                   y[LctState::template get_first_index<InfectionState::InfectedSymptoms>() + i];
            dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>() + i] =
                dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>() + i] - flow;
            dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>() + i + 1] = flow;
        }

        // H'
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] =
            dydt[LctState::template get_first_index<InfectionState::Recovered>()] +
            dydt[LctState::template get_first_index<InfectionState::InfectedSevere>()] *
                (1 - params.template get<SeverePerInfectedSymptoms>());
        dydt[LctState::template get_first_index<InfectionState::InfectedSevere>()] =
            dydt[LctState::template get_first_index<InfectionState::InfectedSevere>()] *
            params.template get<SeverePerInfectedSymptoms>();
        for (size_t i = 0; i < LctState::template get_num_subcompartments<InfectionState::InfectedSevere>(); i++) {
            flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedSevere>() *
                   (1 / params.template get<TimeInfectedSevere>()) *
                   y[LctState::template get_first_index<InfectionState::InfectedSevere>() + i];
            dydt[LctState::template get_first_index<InfectionState::InfectedSevere>() + i] =
                dydt[LctState::template get_first_index<InfectionState::InfectedSevere>() + i] - flow;
            dydt[LctState::template get_first_index<InfectionState::InfectedSevere>() + i + 1] = flow;
        }

        // U'
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] =
            dydt[LctState::template get_first_index<InfectionState::Recovered>()] +
            dydt[LctState::template get_first_index<InfectionState::InfectedCritical>()] *
                (1 - params.template get<CriticalPerSevere>());
        dydt[LctState::template get_first_index<InfectionState::InfectedCritical>()] =
            dydt[LctState::template get_first_index<InfectionState::InfectedCritical>()] *
            params.template get<CriticalPerSevere>();
        for (size_t i = 0; i < LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() - 1;
             i++) {
            flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() *
                   (1 / params.template get<TimeInfectedCritical>()) *
                   y[LctState::template get_first_index<InfectionState::InfectedCritical>() + i];
            dydt[LctState::template get_first_index<InfectionState::InfectedCritical>() + i] =
                dydt[LctState::template get_first_index<InfectionState::InfectedCritical>() + i] - flow;
            dydt[LctState::template get_first_index<InfectionState::InfectedCritical>() + i + 1] = flow;
        }
        // Last flow from U has to be divided between R and D.
        // Must be calculated separately in order not to overwrite the already calculated values ​​for R.
        flow = (ScalarType)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() *
               (1 / params.template get<TimeInfectedCritical>()) *
               y[LctState::template get_first_index<InfectionState::Recovered>() - 1];
        dydt[LctState::template get_first_index<InfectionState::Recovered>() - 1] =
            dydt[LctState::template get_first_index<InfectionState::Recovered>() - 1] - flow;
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] =
            dydt[LctState::template get_first_index<InfectionState::Recovered>()] +
            (1 - params.template get<DeathsPerCritical>()) * flow;
        dydt[LctState::template get_first_index<InfectionState::Dead>()] =
            params.template get<DeathsPerCritical>() * flow;
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
        TimeSeries<ScalarType> compartments_ts((Eigen::Index)InfectionState::Count);
        if (!(LctState::Count == subcompartments_ts.get_num_elements())) {
            log_error("Result does not match infectionState of the Model.");
            Eigen::VectorXd wrong_size = Eigen::VectorXd::Constant((Eigen::Index)InfectionState::Count, -1);
            compartments_ts.add_time_point(-1, wrong_size);
            return compartments_ts;
        }
        Eigen::VectorXd compartments((Eigen::Index)InfectionState::Count);
        for (Eigen::Index i = 0; i < subcompartments_ts.get_num_time_points(); ++i) {
            // Use segment of vector of the result with subcompartments of InfectionState with index j and sum up values of subcompartments.
            compartments[(Eigen::Index)InfectionState::Susceptible] = subcompartments_ts[i][0];
            compartments[(Eigen::Index)InfectionState::Exposed] =
                subcompartments_ts[i]
                    .segment(LctState::template get_first_index<InfectionState::Exposed>(),
                             LctState::template get_num_subcompartments<InfectionState::Exposed>())
                    .sum();
            compartments[(Eigen::Index)InfectionState::InfectedNoSymptoms] =
                subcompartments_ts[i]
                    .segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                             LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
                    .sum();
            compartments[(Eigen::Index)InfectionState::InfectedSymptoms] =
                subcompartments_ts[i]
                    .segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                             LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
                    .sum();
            compartments[(Eigen::Index)InfectionState::InfectedSevere] =
                subcompartments_ts[i]
                    .segment(LctState::template get_first_index<InfectionState::InfectedSevere>(),
                             LctState::template get_num_subcompartments<InfectionState::InfectedSevere>())
                    .sum();
            compartments[(Eigen::Index)InfectionState::InfectedCritical] =
                subcompartments_ts[i]
                    .segment(LctState::template get_first_index<InfectionState::InfectedCritical>(),
                             LctState::template get_num_subcompartments<InfectionState::InfectedCritical>())
                    .sum();
            compartments[(Eigen::Index)InfectionState::Recovered] =
                subcompartments_ts[i][LctState::template get_first_index<InfectionState::Recovered>()];
            compartments[(Eigen::Index)InfectionState::Dead] =
                subcompartments_ts[i][LctState::template get_first_index<InfectionState::Dead>()];

            compartments_ts.add_time_point(subcompartments_ts.get_time(i), compartments);
        }

        return compartments_ts;
    }
};

} // namespace lsecir
} // namespace mio

#endif // LCTSECIR_MODEL_H