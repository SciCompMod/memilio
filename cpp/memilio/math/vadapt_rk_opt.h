/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin J. Kuehn, Daniel Abele, René Schmieding
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
#ifndef VADAPT_RK_OPT_H_
#define VADAPT_RK_OPT_H_

#include "memilio/math/adapt_rk.h"

namespace mio
{

class VRKOptIntegratorCore : public RKIntegratorCore {

private:
    mutable Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> kt_values;
    mutable Eigen::VectorXd yt_eval, ytp1_low;//, ytp1_high; // lower and higher order approximation
    mutable Eigen::ArrayXd eps, error_estimate; // tolerance and estimate used for time step adaption

public:
    VRKOptIntegratorCore() : RKIntegratorCore() {}
    VRKOptIntegratorCore(const double abs_tol, const double rel_tol, const double dt_min, const double dt_max) :
        RKIntegratorCore(abs_tol, rel_tol, dt_min, dt_max)
    {}

    bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                                Eigen::Ref<Eigen::VectorXd> ytp1) const
    {
        double t_eval; // shifted time for evaluating yt
        double dt_new; // updated dt
        
        bool converged = false; // carry for convergence criterion
        bool failed_step_size_adapt = false;

        if (yt_eval.size() == 0) {
            yt_eval.resize(yt.size());
            kt_values.resize(yt.size(), m_tab_final.entries_low.size());
        }

        yt_eval = yt;

        while (!converged && !failed_step_size_adapt) {
            // compute first column of kt, i.e. kt_0 for each y in yt_eval
            f(yt_eval, t, kt_values.col(0));

            for (Eigen::Index i = 1; i < kt_values.cols(); i++) {
                // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
                t_eval = t;
                t_eval += m_tab.entries[i - 1][0] *
                            dt; // t_eval = t + c_i * h // note: line zero of Butcher tableau not stored in array !
                // use ytp1_low as temporary storage for evaluating kt_values[i]
                ytp1_low = yt_eval;
                for (Eigen::VectorXd::Index k = 1; k < m_tab.entries[i - 1].size(); k++) {
                    ytp1_low += (dt * m_tab.entries[i - 1][k]) * kt_values.col(k-1);
                }
                // get the derivatives, i.e., compute kt_i for all y in ytp1_low: kt_i = f(t_eval, ytp1_low)
                f(ytp1_low, t_eval, kt_values.col(i));

            }
            // calculate low order estimate
            ytp1_low = yt_eval;
            ytp1_low += (dt * (kt_values * m_tab_final.entries_low)); // TODO ?combine with assignment?
            // truncation error estimate: yt_low - yt_high = O(h^(p+1)) where p = order of convergence 
            error_estimate = dt * (kt_values * (m_tab_final.entries_high - m_tab_final.entries_low)).array().abs();

            // calculate mixed tolerance
            eps = m_abs_tol + ytp1_low.array().abs() * m_rel_tol;//1e-6;

            converged = (error_estimate <= eps).all(); // convergence criterion
            
            if (converged) {
                // if sufficiently exact, take lower order approximation (higher order is not always higher accuracy!)
                ytp1 = ytp1_low;
                t += dt; // this is the t where ytp1 belongs to
            }
            // else: repeat the calculation (with updated dt)

            // compute new value for dt
            // converged implies eps/error_estimate >= 1, so dt will be increased for the next step
            // hence !converged implies 0 < eps/error_estimate < 1, strictly decreasing dt
            dt_new = dt * std::pow((eps / error_estimate).minCoeff(), (1. / (m_tab_final.entries_low.size() - 1)));
            // safety factor for more conservative step increases,
            // and to avoid dt_new -> dt for step decreases when |error_estimate - eps| -> 0
            dt_new *= 0.9;
            // check if updated dt stays within desired bounds
            // this check should be fine even if dt_new == inf; comparing with NaN should also work, since its always false
            if (m_dt_min < dt_new && dt_new < m_dt_max) {
                dt = dt_new; // update dt for next step
            }
            else {
                failed_step_size_adapt = true;
            }
        }
        return !failed_step_size_adapt;
    }
};

} // namespace mio

#endif