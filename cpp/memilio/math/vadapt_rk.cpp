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
#include "memilio/math/vadapt_rk.h"

namespace mio
{

bool VRKIntegratorCore::step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                             Eigen::Ref<Eigen::VectorXd> ytp1) const
{
    double t_eval; // shifted time for evaluating yt

    bool failed_step_size_adapt = false;
    bool converged              = false;

    if (m_kt_values.size() == 0) {
        m_ytp1_low.resize(yt.size());
        m_ytp1_high.resize(yt.size());
        m_kt_values.resize(yt.size(), m_tab_final.entries_low.size());
    }

    dt *= 2;

    while (!converged && !failed_step_size_adapt) {
        dt *= 0.5;
        // compute first column of kt, i.e. kt_0 for each y in yt
        f(yt, t, m_kt_values.col(0));

        for (Eigen::Index i = 1; i < m_kt_values.cols(); i++) {
            // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
            t_eval = t;
            t_eval += m_tab.entries[i - 1][0] *
                      dt; // t_eval = t + c_i * h // note: line zero of Butcher tableau not stored in array !
            // use m_ytp1_low as temporary storage for evaluating m_kt_values[i]
            m_ytp1_low = yt;
            for (Eigen::VectorXd::Index k = 1; k < m_tab.entries[i - 1].size(); k++) {
                m_ytp1_low += (dt * m_tab.entries[i - 1][k]) * m_kt_values.col(k - 1);
            }
            // get the derivatives, i.e., compute kt_i for all y in m_ytp1_low: kt_i = f(t_eval, m_ytp1_low)
            f(m_ytp1_low, t_eval, m_kt_values.col(i));
        }
        // calculate low order estimate
        m_ytp1_low           = yt + (dt * (m_kt_values * m_tab_final.entries_low));
        m_ytp1_high = yt + (dt * (m_kt_values * m_tab_final.entries_high));

        auto err  = (m_ytp1_low - m_ytp1_high).array().abs();
        auto eps  = dt * (m_abs_tol + m_ytp1_low.array().abs() * m_rel_tol);
        converged = (err <= eps).all();

        if (converged || dt < 2 * m_dt_min + 1e-6) {
            // if sufficiently exact, return m_ytp1_low as lower  order approximation (higher order is not always higher accuracy)
            ytp1 = m_ytp1_low;

            if (dt < 2 * m_dt_min + 1e-6) {
                failed_step_size_adapt = true;
            }
            t += dt; // this is the t where ytp1 belongs to

            const double h = std::pow(0.5, m_tab_final.entries_low.size() - 1);
            if ((err <= h * eps).all() &&
                2 * dt < m_dt_max) { // error of doubled step size is about 2^(p+1) times as large
                dt = 2 * dt;
            }
        }
    }

    return !failed_step_size_adapt;
}

} // namespace mio
