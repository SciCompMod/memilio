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

#ifndef IDE_SEIR_H
#define IDE_SEIR_H

#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "ide_seir/parameters.h"
#include "ide_seir/infection_state.h"
#include "memilio/utils/logging.h"

#include <iostream>

#include <vector>

namespace mio
{
namespace iseir
{
class Model
{
    using Pa  = ParametersBase;
    using Vec = Eigen::VectorX<ScalarType>;

public:
    /**
        * @brief Create an IDE SEIR model.
        *
        * @param[in, out] init TimeSeries with the initial values of the number of susceptibles at associated initial times.
        *   The time steps in this vector should be equidistant and equal to the time step used for the simulation.
        *   A certain history of time steps and values for susceptibles is needed.
        *   A warning is displayed if the condition is violated.
        *   Co be more precise, the first time point needs to be smaller than -(k-1)*TimeStep with
        *       k=ceil((InfectiousTime + LatencyTime)/TimeStep).
        *   The last time point in this vector should be a time 0.
        * @param[in] dt_init The size of the time step used for numerical simulation.
        * @param[in] N_init The population of the considered region.
        */
    Model(TimeSeries<ScalarType>&& init, ScalarType dt_init, int N_init, const Pa& Parameterset_init = Pa())
        : parameters{Parameterset_init}
        , m_result{std::move(init)}
        , m_result_SEIR{TimeSeries<ScalarType>(4)}
        , m_dt{dt_init}
        , m_N{N_init}
    {
    }

    /**
        * @brief Simulate the evolution of infection numbers with the given IDE SEIR model.
        *
        * The simulation is performed by solving the underlying model equation numerically.
        * Here, an integro-differential equation is to be solved. The model parameters and the initial data are used.
        *
        * @param[in] t_max Last simulation day.
        *   If the last point of time of the initial TimeSeries was 0, the simulation will be executed for t_max days.
        * @return The result of the simulation, stored in a TimeSeries with simulation time and
        *       associated number of susceptibles.
        */
    TimeSeries<ScalarType> const& simulate(int t_max)
    {
        m_l = (int)std::floor(parameters.template get<LatencyTime>() / m_dt);
        m_k =
            (int)std::ceil((parameters.template get<InfectiousTime>() + parameters.template get<LatencyTime>()) / m_dt);
        if (m_result.get_time(0) > -(m_k - 1) * m_dt) {
            log_warning("Constraint check: Initial data starts later than necessary. The simulation may be distorted. "
                        "Start the data at time {:.4f} at the latest.",
                        -(m_k - 1) * m_dt);
        }

        while (m_result.get_last_time() < t_max) {
            m_result.add_time_point(m_result.get_last_time() + m_dt);
            Eigen::Index idx = m_result.get_num_time_points();

            // R0t is the effective reproduction number at time t
            auto R0t1 = parameters.template get<ContactFrequency>().get_cont_freq_mat().get_matrix_at(
                            SimulationTime<ScalarType>(m_result.get_time(idx - 2)))(0, 0) *
                        parameters.template get<TransmissionRisk>() * parameters.template get<InfectiousTime>();
            auto R0t2 = parameters.template get<ContactFrequency>().get_cont_freq_mat().get_matrix_at(
                            SimulationTime<ScalarType>(m_result.get_last_time()))(0, 0) *
                        parameters.template get<TransmissionRisk>() * parameters.template get<InfectiousTime>();

            m_result.get_last_value() =
                Vec::Constant(1, m_result[idx - 2][0] * exp(m_dt * (0.5 * 1 / m_N) *
                                                            (R0t1 * num_integration_inner_integral(idx - 2) +
                                                             R0t2 * num_integration_inner_integral(idx - 1))));
        }
        return m_result;
    }

    /**
        * @brief Calculate the distribution of the population in E, I and, R based on the calculated values for S.
        *
        * The values are calculated using the average latency and infection time, not using model equations.
        * The simulated values of S are used for this purpose, so the simulate() function should be called beforehand.
        *
        * @return The result of the calculation stored in an TimeSeries. The TimeSeries contains the simulation time and an
        *   associated Vector with values for S, E, I, and R.
        */
    TimeSeries<ScalarType> const& calculate_EIR()
    {
        Eigen::Index num_points = m_result.get_num_time_points();
        ScalarType S, E, I, R;
        for (Eigen::Index i = m_k; i < num_points; ++i) {
            S = m_result[i][Eigen::Index(InfectionState::S)];
            E = m_result[i - m_l][Eigen::Index(InfectionState::S)] - S;
            I = m_result[i - m_k][Eigen::Index(InfectionState::S)] - m_result[i - m_l][Eigen::Index(InfectionState::S)];
            R = m_N - S - E - I;
            m_result_SEIR.add_time_point(m_result.get_time(i), (Vec(4) << S, E, I, R).finished());
        }
        return m_result_SEIR;
    }

    // Used Parameters for the simulation.
    Pa parameters{};

private:
    /**
        * @brief Density of the generalized beta distribution used for the function f_{beta} of the IDE SEIR model.
        *
        * @param[in] tau evaluation point of the generalized Beta distribution.
        * @param[in] p parameter p of the generalized Beta distribution.
        * @param[in] q parameter q of the generalized Beta distribution.
        * @result Evaluation of the generalized beta distribution at the given evaluation point.
        */
    ScalarType generalized_beta_distribution(ScalarType tau, ScalarType p = 3.0, ScalarType q = 10.0) const
    {
        if ((parameters.template get<LatencyTime>() < tau) &&
            (parameters.template get<InfectiousTime>() + parameters.template get<LatencyTime>() > tau)) {
            return std::tgamma(p + q) * pow(tau - parameters.template get<LatencyTime>(), p - 1) *
                   pow(parameters.template get<InfectiousTime>() + parameters.template get<LatencyTime>() - tau,
                       q - 1) /
                   (std::tgamma(p) * std::tgamma(q) * std::pow(parameters.template get<InfectiousTime>(), p + q - 1));
        }
        return 0.0;
    }

    /**
        * @brief Numerical differentiation of one compartment using a central difference quotient.
        *
        * @param[in] ts_ide TimeSeries with the time steps already calculated.
        *       Used as function values in numerical differentiation.
        * @param[in] idx Time index at which the numerical differentiation should be performed.
        * @param[in] compartment Compartment for which the numerical differentiation is to be performed.
        * @return Numerically approximated derivative of the function belonging to the compartment at the point t[idx].
        */
    ScalarType central_difference_quotient(TimeSeries<ScalarType> const& ts_ide, InfectionState compartment,
                                           Eigen::Index idx) const
    {
        return (ts_ide[idx + 1][Eigen::Index(compartment)] - ts_ide[idx - 1][Eigen::Index(compartment)]) / (2 * m_dt);
    }

    /**
        * @brief Numerical integration of the inner integral of the integro-differential equation for the group S using
        *    a trapezoidal sum.
        *
        * @param[in] idx Index of the point of time used in the inner integral.
        * @return Result of the numerical integration.
        */
    ScalarType num_integration_inner_integral(Eigen::Index idx) const
    {
        ScalarType res = 0.5 * (generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(idx - m_k)) *
                                    central_difference_quotient(m_result, InfectionState::S, m_k) +
                                generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(idx - m_l)) *
                                    central_difference_quotient(m_result, InfectionState::S, idx - m_l));
        Eigen::Index i = idx - m_k + 1;
        while (i <= idx - m_l - 2) {
            res += (generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(i)) *
                    central_difference_quotient(m_result, InfectionState::S, i));
            ++i;
        }
        return res;
    }

    // TimeSeries containing points of time and the corresponding number of susceptibles.
    TimeSeries<ScalarType> m_result;
    // TimeSeries containing points of time and the corresponding number of susceptibles, exposed,
    // infected and recovered.
    TimeSeries<ScalarType> m_result_SEIR;

    // Timestep used for simulation.
    ScalarType m_dt{0};
    // Population of the considered region.
    int m_N{0};

    // Two Indices used for simulation.
    Eigen::Index m_k{0};
    Eigen::Index m_l{0};
};
} // namespace iseir
} // namespace mio
#endif
