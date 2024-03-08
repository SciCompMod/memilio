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
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"
#include "memilio/epidemiology/lct_infection_state.h"

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
template <unsigned int NumExposed, unsigned int NumInfectedNoSymptoms, unsigned int NumInfectedSymptoms,
          unsigned int NumInfectedSevere, unsigned int NumInfectedCritical>
class Model
{

public:
    using LctState =
        LctInfectionState<InfectionState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms, NumInfectedSevere,
                          NumInfectedCritical, 1, 1>; ///< This class specifies the number of subcompartments.

    /**
     * @brief Constructor to create an LCT SECIR Model.
     *
     * @param[in] init Vector with initial values for all infection states inclusive subcompartments.
     * @param[in, out] parameters_init Specifies Parameters necessary for the Model. 
     */
    Model(Eigen::VectorXd init, Parameters&& parameters_init = Parameters())
        : parameters{parameters_init}
        , m_initial_values{std::move(init)}
    {
        m_N0 = m_initial_values.sum();
    }

    /**
     * @brief Checks constraints of the model inclusive check for model parameters.
     */
    bool check_constraints() const
    {
        if (!(LctState::Count == m_initial_values.size())) {
            log_error("Size of the initial values does not match subcompartments.");
            return true;
        }
        for (unsigned int i = 0; i < LctState::Count; i++) {
            if (m_initial_values[i] < 0) {
                log_warning(
                    "Initial values for one subcompartment are less than zero. Simulation results are not realistic.");
                return true;
            }
        }

        return parameters.check_constraints();
    }

    /**
     * @brief Evaulates the right-hand-side f of the LCT dydt = f(y, t).
     *
     * The LCT-SECIR model is defined through ordinary differential equations of the form dydt = f(y, t). 
     * y is a vector containing number of individuals for each (sub-) compartment.
     * This function evaluates the right-hand-side f of the ODE and can be used in an ODE solver.
     * @param[in] y the current state of the model
     * @param[in] t the current time
     * @param[out] dydt a reference to the calculated output
     */
    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorXd> y, ScalarType t, Eigen::Ref<Eigen::VectorXd> dydt) const
    {
        dydt.setZero();

        ScalarType C     = 0;
        ScalarType I     = 0;
        ScalarType dummy = 0;

        // Calculate sum of all subcompartments for InfectedNoSymptoms.
        C = y.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                      LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
                .sum();
        // Calculate sum of all subcompartments for InfectedSymptoms.
        I = y.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                      LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
                .sum();

        // S'
        ScalarType season_val =
            1 + parameters.get<Seasonality>() *
                    sin(3.141592653589793 * (std::fmod((parameters.get<StartDay>() + t), 365.0) / 182.5 + 0.5));
        dydt[0] = -y[0] / (m_N0 - y[LctState::template get_first_index<InfectionState::Dead>()]) * season_val *
                  parameters.get<TransmissionProbabilityOnContact>() *
                  parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(0, 0) *
                  (parameters.get<RelativeTransmissionNoSymptoms>() * C +
                   parameters.get<RiskOfInfectionFromSymptomatic>() * I);

        // E'
        dydt[1] = -dydt[0];
        for (Eigen::Index i = 0; i < LctState::template get_num_subcompartments<InfectionState::Exposed>(); i++) {
            // Dummy stores the value of the flow from dydt[1 + i] to dydt[2 + i].
            // 1+i is always the index of a (sub-)compartment of E and 2+i can also be the index of the first (sub-)compartment of C.
            dummy = LctState::template get_num_subcompartments<InfectionState::Exposed>() *
                    (1 / parameters.get<TimeExposed>()) * y[1 + i];
            // Subtract flow from dydt[1 + i] and add to dydt[2 + i].
            dydt[1 + i] = dydt[1 + i] - dummy;
            dydt[2 + i] = dummy;
        }

        // C'
        for (Eigen::Index i = 0; i < LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>();
             i++) {
            dummy = LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>() *
                    (1 / parameters.get<TimeInfectedNoSymptoms>()) *
                    y[LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() + i];
            dydt[LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() + i] =
                dydt[LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() + i] - dummy;
            dydt[LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() + i + 1] = dummy;
        }

        // I'
        // Flow from last (sub-) compartment of C must be split between I_1 and R.
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] =
            dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>()] *
            parameters.get<RecoveredPerInfectedNoSymptoms>();
        dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>()] =
            dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>()] *
            (1 - parameters.get<RecoveredPerInfectedNoSymptoms>());

        for (Eigen::Index i = 0; i < LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>();
             i++) {
            dummy = LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>() *
                    (1 / parameters.get<TimeInfectedSymptoms>()) *
                    y[LctState::template get_first_index<InfectionState::InfectedSymptoms>() + i];
            dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>() + i] =
                dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>() + i] - dummy;
            dydt[LctState::template get_first_index<InfectionState::InfectedSymptoms>() + i + 1] = dummy;
        }

        // H'
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] =
            dydt[LctState::template get_first_index<InfectionState::Recovered>()] +
            dydt[LctState::template get_first_index<InfectionState::InfectedSevere>()] *
                (1 - parameters.get<SeverePerInfectedSymptoms>());
        dydt[LctState::template get_first_index<InfectionState::InfectedSevere>()] =
            dydt[LctState::template get_first_index<InfectionState::InfectedSevere>()] *
            parameters.get<SeverePerInfectedSymptoms>();
        for (Eigen::Index i = 0; i < LctState::template get_num_subcompartments<InfectionState::InfectedSevere>();
             i++) {
            dummy = LctState::template get_num_subcompartments<InfectionState::InfectedSevere>() *
                    (1 / parameters.get<TimeInfectedSevere>()) *
                    y[LctState::template get_first_index<InfectionState::InfectedSevere>() + i];
            dydt[LctState::template get_first_index<InfectionState::InfectedSevere>() + i] =
                dydt[LctState::template get_first_index<InfectionState::InfectedSevere>() + i] - dummy;
            dydt[LctState::template get_first_index<InfectionState::InfectedSevere>() + i + 1] = dummy;
        }

        // U'
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] =
            dydt[LctState::template get_first_index<InfectionState::Recovered>()] +
            dydt[LctState::template get_first_index<InfectionState::InfectedCritical>()] *
                (1 - parameters.get<CriticalPerSevere>());
        dydt[LctState::template get_first_index<InfectionState::InfectedCritical>()] =
            dydt[LctState::template get_first_index<InfectionState::InfectedCritical>()] *
            parameters.get<CriticalPerSevere>();
        for (Eigen::Index i = 0; i < LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() - 1;
             i++) {
            dummy = LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() *
                    (1 / parameters.get<TimeInfectedCritical>()) *
                    y[LctState::template get_first_index<InfectionState::InfectedCritical>() + i];
            dydt[LctState::template get_first_index<InfectionState::InfectedCritical>() + i] =
                dydt[LctState::template get_first_index<InfectionState::InfectedCritical>() + i] - dummy;
            dydt[LctState::template get_first_index<InfectionState::InfectedCritical>() + i + 1] = dummy;
        }
        // Last flow from U has to be divided between R and D.
        // Must be calculated separately in order not to overwrite the already calculated values ​​for R.
        dummy = LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() *
                (1 / parameters.get<TimeInfectedCritical>()) *
                y[LctState::template get_first_index<InfectionState::Recovered>() - 1];
        dydt[LctState::template get_first_index<InfectionState::Recovered>() - 1] =
            dydt[LctState::template get_first_index<InfectionState::Recovered>() - 1] - dummy;
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] =
            dydt[LctState::template get_first_index<InfectionState::Recovered>()] +
            (1 - parameters.get<DeathsPerCritical>()) * dummy;
        dydt[LctState::template get_first_index<InfectionState::Dead>()] = parameters.get<DeathsPerCritical>() * dummy;
    }

    /**
     * @brief Cumulates a simulation result with subcompartments to produce a result that divides the population only into the infection states defined in InfectionState.
     *
     * If the model is used for simulation, we will get a result in form of a TimeSeries with infection states divided in subcompartments.
     * The function calculates a TimeSeries without subcompartmens from another TimeSeries with subcompartments. 
     * This is done by summing up the numbers in the subcompartments.
     * @param[in] result result of a simulation with the model.
     * @return result of the simulation divided in the Base infection states. 
     *  Returns TimeSeries with values -1 if calculation is not possible.
     */
    TimeSeries<ScalarType> calculate_populations(const TimeSeries<ScalarType>& result) const
    {
        if (!(LctState::Count == result.get_num_elements())) {
            log_error("Result does not match infectionState of the Model.");
            TimeSeries<ScalarType> populations((int)InfectionState::Count);
            Eigen::VectorXd wrong_size = Eigen::VectorXd::Constant((int)InfectionState::Count, -1);
            populations.add_time_point(-1, wrong_size);
            return populations;
        }
        TimeSeries<ScalarType> populations((int)InfectionState::Count);
        Eigen::VectorXd dummy((int)InfectionState::Count);
        for (Eigen::Index i = 0; i < result.get_num_time_points(); ++i) {
            // Use segment of vector of the result with subcompartments of InfectionState with index j and sum up values of subcompartments.
            dummy[(int)InfectionState::Susceptible] = result[i][0];
            dummy[(int)InfectionState::Exposed] =
                result[i]
                    .segment(LctState::template get_first_index<InfectionState::Exposed>(),
                             LctState::template get_num_subcompartments<InfectionState::Exposed>())
                    .sum();
            dummy[(int)InfectionState::InfectedNoSymptoms] =
                result[i]
                    .segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                             LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
                    .sum();
            dummy[(int)InfectionState::InfectedSymptoms] =
                result[i]
                    .segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                             LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
                    .sum();
            dummy[(int)InfectionState::InfectedSevere] =
                result[i]
                    .segment(LctState::template get_first_index<InfectionState::InfectedSevere>(),
                             LctState::template get_num_subcompartments<InfectionState::InfectedSevere>())
                    .sum();
            dummy[(int)InfectionState::InfectedCritical] =
                result[i]
                    .segment(LctState::template get_first_index<InfectionState::InfectedCritical>(),
                             LctState::template get_num_subcompartments<InfectionState::InfectedCritical>())
                    .sum();
            dummy[(int)InfectionState::Recovered] =
                result[i][LctState::template get_first_index<InfectionState::Recovered>()];
            dummy[(int)InfectionState::Dead] = result[i][LctState::template get_first_index<InfectionState::Dead>()];

            populations.add_time_point(result.get_time(i), dummy);
        }

        return populations;
    }

    /**
     * @brief Returns the initial values for the model.
     *
     * This can be used as initial conditions in an ODE solver.
     * @return Vector with initial values for all (sub-)compartments.
     */
    Eigen::VectorXd get_initial_values()
    {
        return m_initial_values;
    }

    /**
     * @brief Sets the initial values for the model.
     *
     * @param[in] init Vector with initial values for all infection states inclusive subcompartments.
     */
    void set_initial_values(Eigen::VectorXd init)
    {
        m_initial_values = init;
        m_N0 = m_initial_values.sum();
    }

    Parameters parameters{}; ///< Parameters of the model.

private:
    Eigen::VectorXd m_initial_values; ///< Initial values of the model.
    ScalarType m_N0{
        0}; ///< Total population size at time t_0 for the considered region (inclusive initial value for Dead).
};

} // namespace lsecir
} // namespace mio

#endif // LCTSECIR_MODEL_H