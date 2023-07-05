/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include <string>

namespace mio
{
namespace lsecir
{
class Model
{
    using ParameterSet = Parameters;

public:
    /**
     * @brief Constructor to create an LCT SECIR model.
     *
     * @param[in] init Vector with initial values for all infection states inclusive subcompartments.
     * @param[in,out] InfectionState_init InfectionStates for the model, specifies number of Subcompartments for each infection state.
     * @param[in, out] Parameterset_init Specifies Parameters necessary for the model. 
     */
    Model(Eigen::VectorXd init, const InfectionState InfectionState_init = InfectionState(),
          const ParameterSet& Parameterset_init = ParameterSet());

    /**
     * @brief Checks constraints of the model inclusive check for model parameters.
     */
    void check_constraints() const;

    /**
     * @brief Evaulates the right-hand-side f of the LCT dydt = f(y, t).
     *
     * The LCT-SECIR model is defined through ordinary differetial equations of the form dydt = f(y, t). 
     * y is a vector containing number of individuals for each (sub-) compartment.
     * This function evaluates the right-hand-side f of the ODE and can be used in an ODE solver.
     * @param[in] y the current state of the model
     * @param[in] t the current time
     * @param[out] dydt a reference to the calculated output
     */
    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorXd> y, ScalarType t,
                              Eigen::Ref<Eigen::VectorXd> dydt) const;

    /**
     * @brief Calculates the population divided in states defined in InfectionStateBase out of a result with subcompartments.
     *
     * If the model is used for simulation, we will get a result in form of a TimeSeries with infection states divided in Subcompartments.
     * Function transforms this TimeSeries in another TimeSeries with just the Basic InfectionStates. 
     * This is done by summing up the numbers in the Subcompartments.
     * @param[in] result result of a simulation with the model.
     * @return result of the simulation divided in the Base infection states.
     */
    TimeSeries<ScalarType> calculate_populations(const TimeSeries<ScalarType>& result) const;

    /**
     * @brief Gives a string with the names of base infection states.
     *
     * Can be used as a heading for printing the result without subcompartments.
     * @return heading in form of a string.
     */
    std::string get_heading_CompartmentsBase() const;

    /**
     * @brief Gives a string with the names of infection stateswith subcompartmens.
     *
     * Can be used as a heading for printing the result with subcompartments.
     * @return heading in form of a string.
     */
    std::string get_heading_Subcompartments() const;

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

    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    InfectionState m_InfectionStates; ///< InfectionState specifies number of subcompartments.

private:
    Eigen::VectorXd m_initial_values; ///< Initial values of the model.
    ScalarType m_N0{0}; ///< Total population size of living people at time t_0 for the considered region.
};

} // namespace lsecir
} // namespace mio

#endif // LCTSECIR_MODEL_H
