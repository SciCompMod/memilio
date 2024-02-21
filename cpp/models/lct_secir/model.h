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
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"

namespace mio
{
namespace lsecir
{

/**
 * @brief Class that defines an LCT-SECIR model.
 *
 * @tparam n_Exposed The number of subcompartents used for the Exposed compartment.
 * @tparam n_InfectedNoSymptoms The number of subcompartents used for the InfectedNoSymptoms compartment. 
 * @tparam n_InfectedSymptoms The number of subcompartents used for the InfectedSymptoms compartment.
 * @tparam n_InfectedSevere The number of subcompartents used for the InfectedSevere compartment.
 * @tparam n_InfectedCritical The number of subcompartents used for the InfectedCritical compartment.
 */
template <unsigned int n_Exposed, unsigned int n_InfectedNoSymptoms, unsigned int n_InfectedSymptoms,
          unsigned int n_InfectedSevere, unsigned int n_InfectedCritical>
class Model
{
    using InfState =
        InfectionState<InfectionStateBase, 1, n_Exposed, n_InfectedNoSymptoms, n_InfectedSymptoms, n_InfectedSevere,
                       n_InfectedCritical, 1, 1>; ///< This class specifies the number of subcompartments.

public:
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
        if (!(InfState::get_count() == m_initial_values.size())) {
            log_error("Size of the initial values does not match subcompartments.");
            return true;
        }
        for (unsigned int i = 0; i < InfState::get_count(); i++) {
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
        C = y.segment(InfState::template get_firstindex<InfectionStateBase::InfectedNoSymptoms>(),
                      InfState::template get_number<InfectionStateBase::InfectedNoSymptoms>())
                .sum();
        // Calculate sum of all subcompartments for InfectedSymptoms.
        I = y.segment(InfState::template get_firstindex<InfectionStateBase::InfectedSymptoms>(),
                      InfState::template get_number<InfectionStateBase::InfectedSymptoms>())
                .sum();

        // S'
        ScalarType season_val =
            1 + parameters.get<Seasonality>() *
                    sin(3.141592653589793 * (std::fmod((parameters.get<StartDay>() + t), 365.0) / 182.5 + 0.5));
        dydt[0] = -y[0] / (m_N0 - y[InfState::template get_firstindex<InfectionStateBase::Dead>()]) * season_val *
                  parameters.get<TransmissionProbabilityOnContact>() *
                  parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(0, 0) *
                  (parameters.get<RelativeTransmissionNoSymptoms>() * C +
                   parameters.get<RiskOfInfectionFromSymptomatic>() * I);

        // E'
        dydt[1] = -dydt[0];
        for (Eigen::Index i = 0; i < InfState::template get_number<InfectionStateBase::Exposed>(); i++) {
            // Dummy stores the value of the flow from dydt[1 + i] to dydt[2 + i].
            // 1+i is always the index of a (sub-)compartment of E and 2+i can also be the index of the first (sub-)compartment of C.
            dummy = InfState::template get_number<InfectionStateBase::Exposed>() * (1 / parameters.get<TimeExposed>()) *
                    y[1 + i];
            // Subtract flow from dydt[1 + i] and add to dydt[2 + i].
            dydt[1 + i] = dydt[1 + i] - dummy;
            dydt[2 + i] = dummy;
        }

        // C'
        for (Eigen::Index i = 0; i < InfState::template get_number<InfectionStateBase::InfectedNoSymptoms>(); i++) {
            dummy = InfState::template get_number<InfectionStateBase::InfectedNoSymptoms>() *
                    (1 / parameters.get<TimeInfectedNoSymptoms>()) *
                    y[InfState::template get_firstindex<InfectionStateBase::InfectedNoSymptoms>() + i];
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedNoSymptoms>() + i] =
                dydt[InfState::template get_firstindex<InfectionStateBase::InfectedNoSymptoms>() + i] - dummy;
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedNoSymptoms>() + i + 1] = dummy;
        }

        // I'
        // Flow from last (sub-) compartment of C must be split between I_1 and R.
        dydt[InfState::template get_firstindex<InfectionStateBase::Recovered>()] =
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSymptoms>()] *
            parameters.get<RecoveredPerInfectedNoSymptoms>();
        dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSymptoms>()] =
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSymptoms>()] *
            (1 - parameters.get<RecoveredPerInfectedNoSymptoms>());

        for (Eigen::Index i = 0; i < InfState::template get_number<InfectionStateBase::InfectedSymptoms>(); i++) {
            dummy = InfState::template get_number<InfectionStateBase::InfectedSymptoms>() *
                    (1 / parameters.get<TimeInfectedSymptoms>()) *
                    y[InfState::template get_firstindex<InfectionStateBase::InfectedSymptoms>() + i];
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSymptoms>() + i] =
                dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSymptoms>() + i] - dummy;
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSymptoms>() + i + 1] = dummy;
        }

        // H'
        dydt[InfState::template get_firstindex<InfectionStateBase::Recovered>()] =
            dydt[InfState::template get_firstindex<InfectionStateBase::Recovered>()] +
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSevere>()] *
                (1 - parameters.get<SeverePerInfectedSymptoms>());
        dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSevere>()] =
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSevere>()] *
            parameters.get<SeverePerInfectedSymptoms>();
        for (Eigen::Index i = 0; i < InfState::template get_number<InfectionStateBase::InfectedSevere>(); i++) {
            dummy = InfState::template get_number<InfectionStateBase::InfectedSevere>() *
                    (1 / parameters.get<TimeInfectedSevere>()) *
                    y[InfState::template get_firstindex<InfectionStateBase::InfectedSevere>() + i];
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSevere>() + i] =
                dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSevere>() + i] - dummy;
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedSevere>() + i + 1] = dummy;
        }

        // U'
        dydt[InfState::template get_firstindex<InfectionStateBase::Recovered>()] =
            dydt[InfState::template get_firstindex<InfectionStateBase::Recovered>()] +
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedCritical>()] *
                (1 - parameters.get<CriticalPerSevere>());
        dydt[InfState::template get_firstindex<InfectionStateBase::InfectedCritical>()] =
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedCritical>()] *
            parameters.get<CriticalPerSevere>();
        for (Eigen::Index i = 0; i < InfState::template get_number<InfectionStateBase::InfectedCritical>() - 1; i++) {
            dummy = InfState::template get_number<InfectionStateBase::InfectedCritical>() *
                    (1 / parameters.get<TimeInfectedCritical>()) *
                    y[InfState::template get_firstindex<InfectionStateBase::InfectedCritical>() + i];
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedCritical>() + i] =
                dydt[InfState::template get_firstindex<InfectionStateBase::InfectedCritical>() + i] - dummy;
            dydt[InfState::template get_firstindex<InfectionStateBase::InfectedCritical>() + i + 1] = dummy;
        }
        // Last flow from U has to be divided between R and D.
        // Must be calculated separately in order not to overwrite the already calculated values ​​for R.
        dummy = InfState::template get_number<InfectionStateBase::InfectedCritical>() *
                (1 / parameters.get<TimeInfectedCritical>()) *
                y[InfState::template get_firstindex<InfectionStateBase::Recovered>() - 1];
        dydt[InfState::template get_firstindex<InfectionStateBase::Recovered>() - 1] =
            dydt[InfState::template get_firstindex<InfectionStateBase::Recovered>() - 1] - dummy;
        dydt[InfState::template get_firstindex<InfectionStateBase::Recovered>()] =
            dydt[InfState::template get_firstindex<InfectionStateBase::Recovered>()] +
            (1 - parameters.get<DeathsPerCritical>()) * dummy;
        dydt[InfState::template get_firstindex<InfectionStateBase::Dead>()] =
            parameters.get<DeathsPerCritical>() * dummy;
    }

    /**
     * @brief Cumulates a simulation result with subcompartments to produce a result that divides the population only into the infection states defined in InfectionStateBase.
     *
     * If the model is used for simulation, we will get a result in form of a TimeSeries with infection states divided in subcompartments.
     * Function transforms this TimeSeries in another TimeSeries with just the Basic infection states. 
     * This is done by summing up the numbers in the subcompartments.
     * @param[in] result result of a simulation with the model.
     * @return result of the simulation divided in the Base infection states. 
     *  Returns TimeSeries with values -1 if calculation is not possible.
     */
    TimeSeries<ScalarType> calculate_populations(const TimeSeries<ScalarType>& result) const
    {
        if (!(InfState::get_count() == result.get_num_elements())) {
            log_error("Result does not match infectionState of the Model.");
            TimeSeries<ScalarType> populations((int)InfectionStateBase::Count);
            Eigen::VectorXd wrong_size = Eigen::VectorXd::Constant((int)InfectionStateBase::Count, -1);
            populations.add_time_point(-1, wrong_size);
            return populations;
        }
        TimeSeries<ScalarType> populations((int)InfectionStateBase::Count);
        Eigen::VectorXd dummy((int)InfectionStateBase::Count);
        for (Eigen::Index i = 0; i < result.get_num_time_points(); ++i) {
            // Use segment of vector of the result with subcompartments of InfectionStateBase with index j and sum up values of subcompartments.
            dummy[(int)InfectionStateBase::Susceptible] = result[i][0];
            dummy[(int)InfectionStateBase::Exposed] =
                result[i]
                    .segment(InfState::template get_firstindex<InfectionStateBase::Exposed>(),
                             InfState::template get_number<InfectionStateBase::Exposed>())
                    .sum();
            dummy[(int)InfectionStateBase::InfectedNoSymptoms] =
                result[i]
                    .segment(InfState::template get_firstindex<InfectionStateBase::InfectedNoSymptoms>(),
                             InfState::template get_number<InfectionStateBase::InfectedNoSymptoms>())
                    .sum();
            dummy[(int)InfectionStateBase::InfectedSymptoms] =
                result[i]
                    .segment(InfState::template get_firstindex<InfectionStateBase::InfectedSymptoms>(),
                             InfState::template get_number<InfectionStateBase::InfectedSymptoms>())
                    .sum();
            dummy[(int)InfectionStateBase::InfectedSevere] =
                result[i]
                    .segment(InfState::template get_firstindex<InfectionStateBase::InfectedSevere>(),
                             InfState::template get_number<InfectionStateBase::InfectedSevere>())
                    .sum();
            dummy[(int)InfectionStateBase::InfectedCritical] =
                result[i]
                    .segment(InfState::template get_firstindex<InfectionStateBase::InfectedCritical>(),
                             InfState::template get_number<InfectionStateBase::InfectedCritical>())
                    .sum();
            dummy[(int)InfectionStateBase::Recovered] =
                result[i][InfState::template get_firstindex<InfectionStateBase::Recovered>()];
            dummy[(int)InfectionStateBase::Dead] =
                result[i][InfState::template get_firstindex<InfectionStateBase::Dead>()];

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

    Parameters parameters{}; ///< Parameters of the model.

private:
    Eigen::VectorXd m_initial_values; ///< Initial values of the model.
    ScalarType m_N0{
        0}; ///< Total population size at time t_0 for the considered region (inclusive initial value for Dead).
};

} // namespace lsecir
} // namespace mio

#endif // LCTSECIR_MODEL_H