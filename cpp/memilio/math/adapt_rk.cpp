/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin J. Kuehn, Daniel Abele, Rene Schmieding
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
#include "memilio/math/adapt_rk.h"

namespace mio
{

TableauFinal::TableauFinal()
{
    entries_low.resize(6);
    entries_high.resize(6);

    entries_low[0]  = 25 / 216.0;
    entries_low[1]  = 0.0;
    entries_low[2]  = 1408 / 2565.0;
    entries_low[3]  = 2197 / 4104.0;
    entries_low[4]  = -0.2;
    entries_low[5]  = 0.0;
    entries_high[0] = 16 / 135.0;
    entries_high[1] = 0.0;
    entries_high[2] = 6656 / 12825.0;
    entries_high[3] = 28561 / 56430.0;
    entries_high[4] = -9 / 50.0;
    entries_high[5] = 2 / 55.0;
}

Tableau::Tableau()
{
    entries.resize(5);
    for (size_t i = 0; i < entries.size(); i++) {
        entries.at(i).resize(i + 2);
    }

    entries[0][0] = 0.25;
    entries[0][1] = 0.25;
    entries[1][0] = 3 / 8.0;
    entries[1][1] = 3 / 32.0;
    entries[1][2] = 9 / 32.0;
    entries[2][0] = 12 / 13.0;
    entries[2][1] = 1932 / 2197.0;
    entries[2][2] = -7200 / 2197.0;
    entries[2][3] = 7296 / 2197.0;
    entries[3][0] = 1.0;
    entries[3][1] = 439 / 216.0;
    entries[3][2] = -8.0;
    entries[3][3] = 3680 / 513.0;
    entries[3][4] = -845 / 4104.0;
    entries[4][0] = 0.5;
    entries[4][1] = -8 / 27.0;
    entries[4][2] = 2.0;
    entries[4][3] = -3544 / 2565.0;
    entries[4][4] = 1859 / 4104.0;
    entries[4][5] = -11 / 40.0;
}

bool RKIntegratorCore::step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                            Eigen::Ref<Eigen::VectorXd> ytp1) const
{
    double t_eval; // shifted time for evaluating yt
    double dt_new; // updated dt

    bool converged              = false; // carry for convergence criterion
    bool failed_step_size_adapt = false;

    if (m_kt_values.size() == 0) {
        m_yt_eval.resize(yt.size());
        m_kt_values.resize(yt.size(), m_tab_final.entries_low.size());
    }

    m_yt_eval = yt;

    while (!converged && !failed_step_size_adapt) {
        // compute first column of kt, i.e. kt_0 for each y in yt_eval
        f(m_yt_eval, t, m_kt_values.col(0));

        for (Eigen::Index i = 1; i < m_kt_values.cols(); i++) {
            // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
            t_eval = t;
            t_eval += m_tab.entries[i - 1][0] *
                      dt; // t_eval = t + c_i * h // note: line zero of Butcher tableau not stored in array
            // use ytp1 as temporary storage for evaluating m_kt_values[i]
            ytp1 = m_yt_eval;
            for (Eigen::VectorXd::Index k = 1; k < m_tab.entries[i - 1].size(); k++) {
                ytp1 += (dt * m_tab.entries[i - 1][k]) * m_kt_values.col(k - 1);
            }
            // get the derivatives, i.e., compute kt_i for all y in ytp1: kt_i = f(t_eval, ytp1_low)
            f(ytp1, t_eval, m_kt_values.col(i));
        }
        // calculate low order estimate
        ytp1 = m_yt_eval;
        ytp1 += (dt * (m_kt_values * m_tab_final.entries_low));
        // truncation error estimate: yt_low - yt_high = O(h^(p+1)) where p = order of convergence
        m_error_estimate = dt * (m_kt_values * (m_tab_final.entries_high - m_tab_final.entries_low)).array().abs();
        // calculate mixed tolerance
        m_eps = m_abs_tol + ytp1.array().abs() * m_rel_tol;

        converged = (m_error_estimate <= m_eps).all(); // convergence criterion

        if (converged) {
            // if sufficiently exact, return ytp1, which currently contains the lower order approximation
            // (higher order is not always higher accuracy)
            t += dt; // this is the t where ytp1 belongs to
        }
        // else: repeat the calculation above (with updated dt)

        // compute new value for dt
        // converged implies eps/error_estimate >= 1, so dt will be increased for the next step
        // hence !converged implies 0 < eps/error_estimate < 1, strictly decreasing dt
        dt_new = dt * std::pow((m_eps / m_error_estimate).minCoeff(), (1. / (m_tab_final.entries_low.size() - 1)));
        // safety factor for more conservative step increases,
        // and to avoid dt_new -> dt for step decreases when |error_estimate - eps| -> 0
        dt_new *= 0.9;
        // check if updated dt stays within desired bounds and update dt for next step
        if (m_dt_min < dt_new) {
            dt = std::min(dt_new, m_dt_max);
        }
        else {
            failed_step_size_adapt = true;
        }
    }
    return !failed_step_size_adapt;
}

} // namespace mio
