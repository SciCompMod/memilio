/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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

#include "ide_seir/model.h"
#include "memilio/utils/logging.h"
#include "ide_seir/parameters.h"
#include "ide_seir/infection_state.h"

#include <iostream>

namespace mio
{
namespace iseir
{
using Pa  = ParametersBase;
using Vec = TimeSeries<double>::Vector;

Model::Model(TimeSeries<double>&& init, double dt_init, int N_init, const Pa& Parameterset_init)
    : parameters{Parameterset_init}
    , m_result{std::move(init)}
    , m_result_SEIR{TimeSeries<double>(4)}
    , m_dt{dt_init}
    , m_N{N_init}
{
}

double Model::generalized_beta_distribution(double tau, double p, double q) const
{
    if ((parameters.get<LatencyTime>() < tau) &&
        (parameters.get<InfectiousTime>() + parameters.get<LatencyTime>() > tau)) {
        return tgamma(p + q) * pow(tau - parameters.get<LatencyTime>(), p - 1) *
               pow(parameters.get<InfectiousTime>() + parameters.get<LatencyTime>() - tau, q - 1) /
               (tgamma(p) * tgamma(q) * pow(parameters.get<InfectiousTime>(), p + q - 1));
    }
    return 0.0;
}

double Model::central_difference_quotient(TimeSeries<double> const& ts_ide, InfectionState compartment,
                                          Eigen::Index idx) const
{
    return (ts_ide[idx + 1][Eigen::Index(compartment)] - ts_ide[idx - 1][Eigen::Index(compartment)]) / (2 * m_dt);
}

double Model::num_integration_inner_integral(Eigen::Index idx) const
{
    double res     = 0.5 * (generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(idx - m_k)) *
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

TimeSeries<double> const& Model::simulate(int t_max)
{

    m_l = (int)std::floor(parameters.get<LatencyTime>() / m_dt);
    m_k = (int)std::ceil((parameters.get<InfectiousTime>() + parameters.get<LatencyTime>()) / m_dt);
    if (m_result.get_time(0) > -(m_k - 1) * m_dt) {
        log_warning("Constraint check: Initial data starts later than necessary. The simulation may be distorted. "
                    "Start the data at time {:.4f} at the latest.",
                    -(m_k - 1) * m_dt);
    }

    while (m_result.get_last_time() < t_max) {
        m_result.add_time_point(m_result.get_last_time() + m_dt);
        Eigen::Index idx = m_result.get_num_time_points();

        // R0t is the effective reproduction number at time t
        auto R0t1 =
            parameters.get<ContactFrequency>().get_cont_freq_mat().get_matrix_at(m_result.get_time(idx - 2))(0, 0) *
            parameters.get<TransmissionRisk>() * parameters.get<InfectiousTime>();
        auto R0t2 =
            parameters.get<ContactFrequency>().get_cont_freq_mat().get_matrix_at(m_result.get_last_time())(0, 0) *
            parameters.get<TransmissionRisk>() * parameters.get<InfectiousTime>();

        m_result.get_last_value() =
            Vec::Constant(1, m_result[idx - 2][0] * exp(m_dt * (0.5 * 1 / m_N) *
                                                        (R0t1 * num_integration_inner_integral(idx - 2) +
                                                         R0t2 * num_integration_inner_integral(idx - 1))));
    }
    return m_result;
}

TimeSeries<double> const& Model::calculate_EIR()
{
    Eigen::Index num_points = m_result.get_num_time_points();
    double S, E, I, R;
    for (Eigen::Index i = m_k; i < num_points; ++i) {
        S = m_result[i][Eigen::Index(InfectionState::S)];
        E = m_result[i - m_l][Eigen::Index(InfectionState::S)] - S;
        I = m_result[i - m_k][Eigen::Index(InfectionState::S)] - m_result[i - m_l][Eigen::Index(InfectionState::S)];
        R = m_N - S - E - I;
        m_result_SEIR.add_time_point(m_result.get_time(i), (Vec(4) << S, E, I, R).finished());
    }
    return m_result_SEIR;
}

void Model::print_result(bool calculated_SEIR) const
{
    if (calculated_SEIR) {
        std::cout << "# time  |  S  |  E  |  I  |  R" << std::endl;
        Eigen::Index num_points = m_result_SEIR.get_num_time_points();
        for (Eigen::Index i = 0; i < num_points; ++i) {
            printf(" %.9f %.9f %.9f %.9f %.9f\n", m_result_SEIR.get_time(i),
                   m_result_SEIR[i][Eigen::Index(InfectionState::S)], m_result_SEIR[i][Eigen::Index(InfectionState::E)],
                   m_result_SEIR[i][Eigen::Index(InfectionState::I)],
                   m_result_SEIR[i][Eigen::Index(InfectionState::R)]);
        }
    }
    else {
        std::cout << "# time  |  number of susceptibles" << std::endl;
        Eigen::Index num_points = m_result.get_num_time_points();
        for (Eigen::Index i = 0; i < num_points; ++i) {
            std::cout << m_result.get_time(i) << "  |  " << m_result[i][Eigen::Index(InfectionState::S)] << std::endl;
        }
    }
}

} // namespace iseir
} // namespace mio
