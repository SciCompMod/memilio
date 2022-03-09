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
#include "memilio/math/adapt_rk_hasty.h"

namespace mio
{

bool HastyRKIntegratorCore::step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                                Eigen::Ref<Eigen::VectorXd> ytp1) const
{
    double conv_crit = 0, max_err = 1, t_eval; // shifted time for evaluating yt

    bool failed_step_size_adapt = false;

    if (m_kt_values.size() == 0) {
        m_yt_eval.resize(yt.size());
        m_kt_values.resize(yt.size(), m_tab_final.entries_low.size());
    }
    m_yt_eval = yt;

    dt *= 2;
    while (max_err > conv_crit && !failed_step_size_adapt) {
        dt *= 0.5;
        // compute first column of kt, i.e. kt_0 for each y in yt
        f(m_yt_eval, t, m_kt_values.col(0));

        for (Eigen::Index i = 1; i < m_kt_values.cols(); i++) {
            // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
            t_eval = t;
            t_eval += m_tab.entries[i - 1][0] *
                      dt; // t_eval = t + c_i * h // note: line zero of Butcher tableau not stored in array !
            // use ytp1 as temporary storage for evaluating m_kt_values[i]
            ytp1 = m_yt_eval;
            for (Eigen::VectorXd::Index k = 1; k < m_tab.entries[i - 1].size(); k++) {
                ytp1 += (dt * m_tab.entries[i - 1][k]) * m_kt_values.col(k - 1);
            }
            // get the derivatives, i.e., compute kt_i for all y in ytp1: kt_i = f(t_eval, ytp1)
            f(ytp1, t_eval, m_kt_values.col(i));
        }
        // calculate low and high order estimate
        ytp1           = m_yt_eval + (dt * (m_kt_values * m_tab_final.entries_low));
        auto ytp1_high = m_yt_eval + (dt * (m_kt_values * m_tab_final.entries_high));

        max_err        = (ytp1 - ytp1_high).lpNorm<Eigen::Infinity>();
        double max_val = ytp1.lpNorm<Eigen::Infinity>();

        conv_crit = dt * (m_abs_tol + max_val * m_rel_tol); // scale conv_crit by dt instead of max_err by 1/dt

        if (max_err <= conv_crit || dt < 2 * m_dt_min + 1e-6) {
            // if sufficiently exact, return ytp1 as lower order approximation (higher order is not always higher accuracy)

            if (dt < 2 * m_dt_min + 1e-6) {
                failed_step_size_adapt = true;
            }
            t += dt; // this is the t where ytp1 belongs to

            const double h = std::pow(0.5, m_tab_final.entries_low.size() - 1);
            if (max_err <= h * conv_crit &&
                2 * dt < m_dt_max) { // error of doubled step size is about 2^(p+1) times as large
                dt = 2 * dt;
            }
        }
    }
    return !failed_step_size_adapt;
}

} // namespace mio
