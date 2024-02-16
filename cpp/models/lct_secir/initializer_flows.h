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

#ifndef LCTSECIR_INITIALIZER_H
#define LCTSECIR_INITIALIZER_H

#include "memilio/config.h"
#include "lct_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "lct_secir/parameters.h"

namespace mio
{
namespace lsecir
{

/**
 * @brief Class that can be used to compute an initialization vector out of flows for a LCT Model.
 */
class Initializer
{
public:
    /**
     * @brief Constructs a new Initializer object.
     *
     * @param[in] flows Initalizing TimeSeries with flows fitting to these defined in InfectionTransition. 
     *      Timesteps should be equidistant.
     * @param[in, out] infectionState_init InfectionState%s for the Initializer, specifies number of Subcompartments for each infection state.
     * @param[in, out] parameters_init Specifies parameters necessary for the Initializer.
     */
    Initializer(TimeSeries<ScalarType>&& flows, InfectionState infectionState_init = InfectionState(),
                Parameters&& parameters_init = Parameters())
        : parameters(parameters_init)
        , infectionState(infectionState_init)
        , m_flows(std::move(flows))
    {
        m_dt = m_flows.get_time(1) - m_flows.get_time(0);
    }

    /**
     * @brief Core function of Initializer.
     *
     * Computes a vector that can be used for the initalization of an LCT model with the number of persons for each subcompartment.
     *
     * @param[in] total_population The total size of the considered population.
     * @param[in] deaths Number of deceased people from the disease at time 0.
     * @param[in] total_confirmed_cases Total number of confirmed cases at time 0.
     * @return Vector with a possible initialization for an LCT model computed out of the flows. 
     *      A subcompartment is set to -1 if calculation was not possible.
     */
    Eigen::VectorXd compute_initializationvector(ScalarType total_population, ScalarType deaths,
                                                 ScalarType total_confirmed_cases) const;

    /**
     * @brief Checks constraints of the Initializer inclusive check for parameters.
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false. 
     */
    bool check_constraints() const;

    /**
     * @brief Setter for the tolerance used to calculate the maximum support of ErlangDensity%s.
     *
     * @param[in] new_tol New tolerance.
     */
    void set_tol_for_support_max(ScalarType new_tol)
    {
        m_tol = new_tol;
    }

    Parameters parameters{}; ///< Parameters with the Initalizers Parameters.
    InfectionState infectionState; ///< InfectionState specifies number of subcompartments.

private:
    /**
     * @brief Computes a vector with the number of people in each compartment for one InfectionStateBase.
     *
     * With this function, partial result of compute_initializationvector are achieved.
     *
     * @param[in] base The InfectionStateBase for which the partial result should be calculated.
     * @param[in] idx_incoming_flow Index of the flow which is relevant for the calculation, so the flow to the InfectionStateBase base.
     * @param[in] transition_rate Specifies the transition rate of the InfectionsStateBase. Is equal to 1 / (expected Time in base).
     * @return Vector with a possible initialization for the subcompartments of base. 
     *      Subcompartment is set to -1 if calculation was not possible.
     */
    Eigen::VectorXd compute_compartment(InfectionStateBase base, Eigen::Index idx_incoming_flow,
                                        ScalarType transition_rate) const;

    TimeSeries<ScalarType> m_flows; ///< TimeSeries with the flows which are used to calculate the initial vector.
    ScalarType m_dt{}; ///< Step size of the times in m_flows and time step for the approximation of the integral.
    ScalarType m_tol{1e-10}; ///< Tolerance used to calculate the maximum support of the ErlangDensity%s.
};

} // namespace lsecir
} // namespace mio

#endif