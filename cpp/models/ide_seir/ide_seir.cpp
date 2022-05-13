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

#include "ide_seir/ide_seir.h"
#include "memilio/utils/logging.h"
#include "ide_seir/parameters.h"

#include <iostream>

/* TODO 
    Test schreiben
    beschreibungen und kommentare rewriten
    */

namespace mio
{
namespace iseir
{
    using Pa  = IdeSeirParameters;
    using Vec = TimeSeries<double>::Vector;

    IdeSeirModel::IdeSeirModel(TimeSeries<double>&& init, double dt_init, int N_init, Pa Parameterset_init)
        : m_parameters{Parameterset_init}
        , m_result{std::move(init)}
        , m_dt{dt_init}
        , m_N{N_init}
    {
    }

    double IdeSeirModel::generalized_beta_distribution(double tau, double p, double q) const
    {
        if ((m_parameters.get<LatencyTime>() < tau) &&
            (m_parameters.get<InfectiousTime>() + m_parameters.get<LatencyTime>() > tau)) {
            return tgamma(p + q) * pow(tau - m_parameters.get<LatencyTime>(), p - 1) *
                   pow(m_parameters.get<InfectiousTime>() + m_parameters.get<LatencyTime>() - tau, q - 1) /
                   (tgamma(p) * tgamma(q) * pow(m_parameters.get<InfectiousTime>(), p + q - 1));
        }
        return 0.0;
    }

    double IdeSeirModel::central_difference_quotient(TimeSeries<double> const& ts_ide, iSeirInfType compartment,
                                                     Eigen::Index idx) const
    {
        return (ts_ide[idx + 1][Eigen::Index(compartment)] - ts_ide[idx - 1][Eigen::Index(compartment)]) / (2 * m_dt);
    }

    double IdeSeirModel::num_integration_inner_integral(Eigen::Index idx) const
    {
        double res     = 0.5 * (generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(idx - m_k)) *
                                central_difference_quotient(m_result, iSeirInfType::S, m_k) +
                            generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(idx - m_l)) *
                                central_difference_quotient(m_result, iSeirInfType::S, idx - m_l));
        Eigen::Index i = idx - m_k + 1;
        while (i <= idx - m_l - 2) {
            res += (generalized_beta_distribution(m_result.get_time(idx) - m_result.get_time(i)) *
                    central_difference_quotient(m_result, iSeirInfType::S, i));
            ++i;
        }
        return res;
    }

    TimeSeries<double> const& IdeSeirModel::simulate(int t_max)
    {

        m_l = (int)std::floor(m_parameters.get<LatencyTime>() / m_dt);
        m_k = (int)std::ceil((m_parameters.get<InfectiousTime>() + m_parameters.get<LatencyTime>()) / m_dt);
        if (m_result.get_time(0) > -(m_k - 1) * m_dt) {
            log_warning("Constraint check: Initial data starts later than necessary. The simulation may be distorted. "
                        "Start the data at time {:.4f} at the latest.",
                        -(m_k - 1) * m_dt);
        }

        while (m_result.get_last_time() < t_max) {
            m_result.add_time_point(m_result.get_last_time() + m_dt);
            Eigen::Index idx = m_result.get_num_time_points();

            // R0t is the effective reproduction number at time t
            auto R0t1 = m_parameters.get<ContactFrequency>().get_cont_freq_mat().get_matrix_at(
                            m_result.get_time(idx - 2))(0, 0) *
                        m_parameters.get<TransmissionRisk>() * m_parameters.get<InfectiousTime>();
            auto R0t2 =
                m_parameters.get<ContactFrequency>().get_cont_freq_mat().get_matrix_at(m_result.get_last_time())(0, 0) *
                m_parameters.get<TransmissionRisk>() * m_parameters.get<InfectiousTime>();

            m_result.get_last_value() =
                Vec::Constant(1, m_result[idx - 2][0] * exp(m_dt * (0.5 * 1 / m_N) *
                                                            (R0t1 * num_integration_inner_integral(idx - 2) +
                                                             R0t2 * num_integration_inner_integral(idx - 1))));
        }
        return m_result;
    }

    TimeSeries<double> const& IdeSeirModel::calculate_EIR()
    {
        Eigen::Index num_points = m_result.get_num_time_points();
        double S, E, I, R;
        for (Eigen::Index i = m_k; i < num_points; ++i) {
            S = m_result[i][Eigen::Index(iSeirInfType::S)];
            E = m_result[i - m_l][Eigen::Index(iSeirInfType::S)] - S;
            I = m_result[i - m_k][Eigen::Index(iSeirInfType::S)] - m_result[i - m_l][Eigen::Index(iSeirInfType::S)];
            R = m_N - S - E - I;
            m_result_SEIR.add_time_point(m_result.get_time(i), (Vec(4) << S, E, I, R).finished());
        }
        return m_result_SEIR;
    }

    void IdeSeirModel::print_result(bool calculated_SEIR) const
    {
        if (calculated_SEIR) {
            std::cout << "time  |  S  |  E  |  I  |  R" << std::endl;
            Eigen::Index num_points = m_result_SEIR.get_num_time_points();
            for (Eigen::Index i = 0; i < num_points; ++i) {
                std::cout << m_result_SEIR.get_time(i) << "  |  " << m_result_SEIR[i][Eigen::Index(iSeirInfType::S)]
                          << "  |  " << m_result_SEIR[i][Eigen::Index(iSeirInfType::E)] << "  |  "
                          << m_result_SEIR[i][Eigen::Index(iSeirInfType::I)] << "  |  "
                          << m_result_SEIR[i][Eigen::Index(iSeirInfType::R)] << std::endl;
            }
        }
        else {
            std::cout << "time  |  number of susceptibles" << std::endl;
            Eigen::Index num_points = m_result.get_num_time_points();
            for (Eigen::Index i = 0; i < num_points; ++i) {
                std::cout << m_result.get_time(i) << "  |  " << m_result[i][Eigen::Index(iSeirInfType::S)] << std::endl;
            }
        }
    }

} // namespace iseir
} // namespace mio