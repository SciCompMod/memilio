/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin J. Kuehn, Daniel Abele
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
#include "epidemiology/math/adapt_rk.h"
#include <epidemiology/utils/logging.h>

namespace epi
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

    double max_err   = 1e10;
    double conv_crit = 1e9;

    std::vector<Eigen::VectorXd> kt_values;

    bool failed_step_size_adapt = false;

    dt        = 2 * dt;
    bool cond = false;

    while (!cond && !failed_step_size_adapt) {
        dt = 0.5 * dt;

        kt_values.resize(0); // remove data from previous loop
        kt_values.resize(m_tab_final.entries_low.size()); // these are the k_ni per y_t, used to compute y_t+1

        for (size_t i = 0; i < kt_values.size();
             i++) { // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
            kt_values[i].resizeLike(yt); // note: yt contains more than one variable since we solve a system of ODEs

            if (i == 0) { // for i==0, we have kt_i=f(t,y(t))
                f(yt, t, kt_values[i]);
                // printf("\n\n t = %.4f\t kn1 = %.4f", t, kt_values[i][0]);
            }
            else {

                double t_eval = t;

                t_eval += m_tab.entries[i - 1][0] *
                          dt; // t_eval = t + c_i * h // note: line zero of Butcher tableau not stored in array !

                // printf("\n t = %.4f", t_eval);

                auto yt_eval = yt.eval();

                // go through the variables of the system: S, E, I, ....
                for (int j = 0; j < yt_eval.size(); j++) {
                    // go through the different kt_1, kt_2, ..., kt_i-1 and add them onto yt: y_eval = yt + h * \sum_{j=1}^{i-1} a_{i,j} kt_j
                    for (size_t k = 1; k < m_tab.entries[i - 1].size(); k++) {
                        yt_eval[j] +=
                            (dt * m_tab.entries[i - 1][k] *
                             kt_values
                                 [k - 1]
                                 [j]); // note the shift in k and k-1 since the first column of 'tab' corresponds to 'b_i' and 'a_ij' starts with the second column
                    }
                }

                // get the derivatives, i.e., compute kt_i for all y at yt_eval: kt_i = f(t_eval, yt_eval)
                f(yt_eval, t_eval, kt_values[i]);

                // printf("\t kn%d = %.4f", i+1, kt_values[i][0]);
            }
        }

        // copy actual yt to compare both new iterates
        auto ytp1_low  = yt.eval(); // lower order approximation (e.g., order 4)
        auto ytp1_high = yt.eval(); // higher order approximation (e.g., order 5)

        std::vector<double> err(yt.size(), 0);
        double max_val = 0;
        max_err        = 0;

        for (int i = 0; i < yt.size(); i++) {

            for (size_t j = 0; j < kt_values.size();
                 j++) { // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
                ytp1_low[i] += (dt * m_tab_final.entries_low[j] * kt_values[j][i]);
                ytp1_high[i] += (dt * m_tab_final.entries_high[j] * kt_values[j][i]);
            }

            err[i] =
                1 / dt *
                std::abs(ytp1_low[i] -
                         ytp1_high[i]); // divide by h=dt since the local error is one order higher than the global one

            // printf("\n i: %lu,\t err_abs %e, err_rel %e\n", i,  err[i], max_val);
            if (err[i] > max_err) {
                max_err = err[i];
            }
            if (max_val < std::abs(ytp1_low[i])) {
                max_val = std::abs(ytp1_low[i]);
            }
        }

        // printf("\n low %.8f\t high %.8f\t max_err: %.8f,\t val: %.8f,\t crit: %.8f,\t %d t: %.8f\t dt: %.8f ",
        //    ytp1_low[0], ytp1_high[0], max_err, max_val, m_abs_tol + max_val * m_rel_tol,
        //    max_err <= (m_abs_tol + max_val * m_rel_tol), t, dt);
        // sleep(1);
        conv_crit = m_abs_tol + max_val * m_rel_tol;

        cond = true;
        for (int i = 0; i < yt.size(); i++) {
            cond = cond && (err[i] <= m_abs_tol + std::abs(ytp1_low[i]) * m_rel_tol);
        }
        if (cond || dt < 2 * m_dt_min + 1e-6) {
            // if sufficiently exact, take 4th order approximation (do not take 5th order : Higher order is not always higher accuracy!)
            ytp1 = ytp1_low;

            for (int i = 0; i < ytp1.size(); i++) {
                if (ytp1[i] < -0.1) {
                    log_warning("Compartment {:d} for Age Group {:d} is negative at time {:0.4f}: {:0.4f}", i % 26,
                                (int)i / 26, t, ytp1[i]);
                }
            }

            if (dt < 2 * m_dt_min + 1e-6) {
                failed_step_size_adapt = true;
            }

            // printf("\n succ: %d t: %.2e dt: %.2e \t ", !failed_step_size_adapt, t, dt);
            // double sum = 0;
            // for (size_t kkk = 0; kkk < (size_t)ytp1.size(); kkk++) {
            //     printf(" [%lu] %.2e ", kkk, ytp1[kkk]);
            //     sum += ytp1[kkk];
            // }
            // printf(" %.2e ", sum);

            t += dt; // this is the t where ytp1 belongs to

            if (max_err <= 0.03 * conv_crit &&
                2 * dt < m_dt_max) { // error of doubled step size is about 32 times as large
                dt = 2 * dt;
            }
        }
    }

    return !failed_step_size_adapt;
}

} // namespace epi
