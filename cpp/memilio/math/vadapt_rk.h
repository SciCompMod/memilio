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
#ifndef VADAPT_RK_H_
#define VADAPT_RK_H_

#include "memilio/math/adapt_rk.h"

namespace mio
{

class VRKIntegratorCore : public RKIntegratorCore {

public:
    VRKIntegratorCore() : RKIntegratorCore() {}
    VRKIntegratorCore(const double abs_tol, const double rel_tol, const double dt_min, const double dt_max) :
        RKIntegratorCore(abs_tol, rel_tol, dt_min, dt_max) {}

    bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                                Eigen::Ref<Eigen::VectorXd> ytp1) const
    {
        double conv_crit = 0, max_err = 1, t_eval; // shifted time for evaluating yt
        double dt_new; // updated dt

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> kt_values(yt.size(), m_tab_final.entries_low.size());

        Eigen::VectorXd ytp1_low, ytp1_high; // lower and higher order approximation (e.g., order 4)

        bool failed_step_size_adapt = false;
        //bool crit2 = false;

        dt *= 2;

        while (max_err > conv_crit && !failed_step_size_adapt) {
        //while (!crit2 && !failed_step_size_adapt) {
            dt *= 0.5;
            // compute first column of kt, i.e. kt_0 for each y in yt
            f(yt, t, kt_values.col(0));

            for (Eigen::Index i = 1; i < kt_values.cols(); i++) {
                // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
                t_eval = t;
                t_eval += m_tab.entries[i - 1][0] *
                            dt; // t_eval = t + c_i * h // note: line zero of Butcher tableau not stored in array !
                // use ytp1_low as temporary storage for evaluating kt_values[i]
                ytp1_low = yt;
                for (Eigen::VectorXd::Index k = 1; k < m_tab.entries[i - 1].size(); k++) {
                    ytp1_low += (dt * m_tab.entries[i - 1][k]) * kt_values.col(k-1);
                }
                // get the derivatives, i.e., compute kt_i for all y in ytp1_low: kt_i = f(t_eval, ytp1_low)
                f(ytp1_low, t_eval, kt_values.col(i));

            }
            // calculate low order estimate
            ytp1_low  = yt + (dt * (kt_values * m_tab_final.entries_low ));
            ytp1_high = yt + (dt * (kt_values * m_tab_final.entries_high));
            
            max_err = (ytp1_low - ytp1_high).lpNorm<Eigen::Infinity>();
            double max_val = ytp1_low.lpNorm<Eigen::Infinity>();

            conv_crit = dt * (m_abs_tol + max_val * m_rel_tol); // scale conv_crit by dt instead of max_err by 1/dt
            //crit2 = ((ytp1_low - ytp1_high).array().abs() <= dt * (m_abs_tol + ytp1_low.array().abs() * m_rel_tol)).all();

            if (max_err <= conv_crit || dt < 2 * m_dt_min + 1e-6) {
            //if (crit2 || dt < 2 * m_dt_min + 1e-6) {
                // if sufficiently exact, take 4th order approximation (do not take 5th order : Higher order is not always higher accuracy!)
                ytp1 = ytp1_low;

                if (dt < 2 * m_dt_min + 1e-6) {
                    failed_step_size_adapt = true;
                }
                t += dt; // this is the t where ytp1 belongs to
                //bool crit3 = ((ytp1_low - ytp1_high).array().abs()  <=  0.03 * dt * (m_abs_tol + ytp1_low.array().abs() * m_rel_tol)).all();

                if (max_err <= 0.03 * conv_crit &&
                //if (crit3 &&
                    2 * dt < m_dt_max) { // error of doubled step size is about 32 times as large
                    dt = 2 * dt;
                }
            }
        }

        return !failed_step_size_adapt;
    }
};

} // namespace mio

#endif