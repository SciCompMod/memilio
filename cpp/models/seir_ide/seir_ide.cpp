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

#include "seir_ide/seir_ide.h"
#include "memilio/epidemiology/contact_matrix.h"
#include <iostream>

namespace mio
{
/* TODO east coast const
    logging in constructor (include)
    constructor umformulieren ( no we need, ->)
    Test schreiben*/

using Vec = TimeSeries<double>::Vector;

IdeModel::IdeModel(TimeSeries<double>&& init, double dt_init, int N_init)
    : m_result{std::move(init)}
    , m_dt{dt_init}
    , m_N{N_init}
{
    m_l = (int)std::floor(m_latency_time / m_dt);
    m_k = (int)std::ceil((m_infectious_time + m_latency_time) / m_dt);
}

void IdeModel::set_latency_time(double latency)
{
    m_latency_time = latency;
    // initialize parameters m_k and m_l used for numerical calculations in the simulation
    m_l = (int)std::floor(m_latency_time / m_dt);
    m_k = (int)std::ceil((m_infectious_time + m_latency_time) / m_dt);
}

void IdeModel::set_infectious_time(double infectious)
{
    m_infectious_time = infectious;
    m_k               = (int)std::ceil((m_infectious_time + m_latency_time) / m_dt);
}

ContactMatrix& IdeModel::get_contact_matrix()
{
    return m_contact_matrix;
}

double IdeModel::generalized_beta_distribution(double tau, double p, double q) const
{
    if ((m_latency_time < tau) && (m_infectious_time + m_latency_time > tau)) {
        return tgamma(p + q) * pow(tau - m_latency_time, p - 1) * pow(m_infectious_time + m_latency_time - tau, q - 1) /
               (tgamma(p) * tgamma(q) * pow(m_infectious_time, p + q - 1));
    }
    return 0.0;
}

double IdeModel::central_difference_quotient(Eigen::Index compartment, Eigen::Index idx) const
{
    return (m_result[idx + 1][compartment] - m_result[idx - 1][compartment]) / (2 * m_dt);
}

double IdeModel::num_integration_inner_integral(Eigen::Index idx) const
{
    double res     = 0.5 * (generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(idx - m_k)) *
                            central_difference_quotient(0, m_k) +
                        generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(idx - m_l)) *
                            central_difference_quotient(0, idx - m_l));
    Eigen::Index i = idx - m_k + 1;
    while (i <= idx - m_l - 2) {
        res += (generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(i)) *
                central_difference_quotient(0, i));
        ++i;
    }
    return res;
}

const TimeSeries<double>& IdeModel::simulate(int t_max)
{
    while (m_result.get_last_time() < t_max) {
        m_result.add_time_point(round((m_result.get_last_time() + m_dt) * 10000.0) / 10000.0);
        Eigen::Index idx = m_result.get_num_time_points();

        //R0t = effective contactfrequency at time t * timeinfectious
        auto R0t1 = m_contact_matrix.get_matrix_at(m_result.get_time(idx - 2))(0, 0) * m_infectious_time;
        auto R0t2 = m_contact_matrix.get_matrix_at(m_result.get_last_time())(0, 0) * m_infectious_time;

        m_result.get_last_value() =
            Vec::Constant(1, m_result[idx - 2][0] * exp(m_dt * (0.5 * 1 / m_N) *
                                                        (R0t1 * num_integration_inner_integral(idx - 2) +
                                                         R0t2 * num_integration_inner_integral(idx - 1))));
    }
    return m_result;
}

const TimeSeries<double>& IdeModel::calculate_EIR()
{
    Eigen::Index num_points = m_result.get_num_time_points();
    double S, E, I, R;
    for (Eigen::Index i = m_k; i < num_points; ++i) {
        S = m_result[i][0];
        E = m_result[i - m_l][0] - S;
        I = m_result[i - m_k][0] - m_result[i - m_l][0];
        R = m_N - S - E - I;
        m_result_SEIR.add_time_point(m_result.get_time(i), (Vec(4) << S, E, I, R).finished());
    }
    return m_result_SEIR;
}

void IdeModel::print_result(bool calculated_SEIR) const
{
    if (calculated_SEIR) {
        std::cout << "time  |  S  |  E  |  I  |  R" << std::endl;
        Eigen::Index num_points = m_result_SEIR.get_num_time_points();
        for (Eigen::Index i = 0; i < num_points; ++i) {
            std::cout << m_result_SEIR.get_time(i) << "  |  " << m_result_SEIR[i][0] << "  |  " << m_result_SEIR[i][1]
                      << "  |  " << m_result_SEIR[i][2] << "  |  " << m_result_SEIR[i][3] << std::endl;
        }
    }
    else {
        std::cout << "time  |  number of susceptibles" << std::endl;
        Eigen::Index num_points = m_result.get_num_time_points();
        for (Eigen::Index i = 0; i < num_points; ++i) {
            std::cout << m_result.get_time(i) << "  |  " << m_result[i][0] << std::endl;
        }
    }
}

} // namespace mio