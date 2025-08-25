/*
* Copyright (C) 2020-2025 MEmilio
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
#ifndef ADAPT_RK_H_
#define ADAPT_RK_H_

#include "memilio/math/integrator.h"

#include <cstdio>
#include <vector>

namespace mio
{

/**
 * Two scheme Runge-Kutta numerical integrator with adaptive step width
 * for ODE y'(t) = f(t,y) which is given by
 *    y_{n+1} = y_n + h*\sum_{i=1}^lb_ik_{ni}
 * with
 *    k_{ni} = f(t_n + c_i*h, y_n + h*\sum_{j=1}^{i-1}a_{ij}k_{nj})
 * where the general Butcher tableau is
 * 0   |
 * c_2 | a_{21}
 * c_3 | a_{31} a_{32}
 * ... | ............
 * c_s | a{s,1} a_{s,2}  ... a_{s,s-1}
 * ----------------------------------------------
 *     | b_1    b_2    ...   b_{s-1}     b_s
 *
 *
 * Actually, only the RK-Fehlberg 4 (5) method is implemented.
 *
 * The RKF4's Butcher tableau is
 * 0     |
 * 1/4   | 1/4
 * 3/8   | 3/32        9/32
 * 12/13 | 1932/2197   −7200/2197    7296/2197
 * 1     | 439/216     -8            3680/513     -845/4104
 * 1/2   | -8/27       2             -3544/2565   1859/4104     -11/40
 * ----------------------------------------------------------------------------
 *       | 25/216      0             1408/2565    2197/4104     -1/5        0
 * The higher order (5th) approximation only differs by the last line which is
 *       | 16/135      0             6656/12825   28561/56430   -9/50       2/55
 *
 */
template <typename FP>
class Tableau
{
public:
    std::vector<Eigen::VectorX<FP>> entries;

    /**
     * @brief default is Runge-Kutta-Fehlberg4(5) tableau
     */
    Tableau()
    {
        entries.resize(5);
        for (size_t i = 0; i < entries.size(); i++) {
            entries[i].resize(i + 2);
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
};

/*
 * Final row(s) of the Butcher tableau
 *
 * Repitition from above: Runge-Kutta-4 is:
 *       | 25/216      0             1408/2565    2197/1404     -1/5        0
 * The higher order (5th) approximation only differs by the last line which is
 *       | 16/135      0             6656/12825   28561/56430   -9/50       2/55
 *
 */
template <typename FP>
class TableauFinal
{
public:
    Eigen::VectorX<FP> entries_low;
    Eigen::VectorX<FP> entries_high;

    /**
     * @brief default is Runge-Kutta-Fehlberg4(5) tableau
     */
    TableauFinal()
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
};

/**
 * @brief Two scheme Runge-Kutta numerical integrator with adaptive step width
 *
 * This class integrates a system of ODEs via the step method
 */
template <typename FP>
class RKIntegratorCore : public OdeIntegratorCore<FP>
{
public:
    /**
     * @brief Setting up the integrator
     */
    RKIntegratorCore()
        : OdeIntegratorCore<FP>(std::numeric_limits<ScalarType>::min(), std::numeric_limits<ScalarType>::max())
        , m_abs_tol(1e-10)
        , m_rel_tol(1e-5)
    {
    }

    /**
     * @brief Set up the integrator
     * @param abs_tol absolute tolerance
     * @param rel_tol relative tolerance
     * @param dt_min lower bound for time step dt
     * @param dt_max upper bound for time step dt
     */
    RKIntegratorCore(const FP abs_tol, const FP rel_tol, const FP dt_min, const FP dt_max)
        : OdeIntegratorCore<FP>(dt_min, dt_max)
        , m_abs_tol(abs_tol)
        , m_rel_tol(rel_tol)
    {
    }

    /// @param tol the required absolute tolerance for the comparison with the Fehlberg approximation
    void set_abs_tolerance(FP tol)
    {
        m_abs_tol = tol;
    }

    /// @param tol the required relative tolerance for the comparison with the Fehlberg approximation
    void set_rel_tolerance(FP tol)
    {
        m_rel_tol = tol;
    }

    /// @param dt_min sets the minimum step size
    void set_dt_min(FP dt_min)
    {
        this->get_dt_min() = dt_min;
    }

    /// @param dt_max sets the maximum step size
    void set_dt_max(FP dt_max)
    {
        this->get_dt_max() = dt_max;
    }

    // Allow setting different RK tableau schemes
    void set_tableaus(const Tableau<FP>& tab, const TableauFinal<FP>& final_tab)
    {
        m_tab       = tab;
        m_tab_final = final_tab;
        m_kt_values.resize(Eigen::NoChange, m_tab_final.entries_low.size());
    }

    /**
     * @brief Make a single integration step of a system of ODEs and adapt the step size
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time
     * @param[in,out] dt current time step size h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
    bool step(const DerivFunction<FP>& f, Eigen::Ref<const Eigen::VectorX<FP>> yt, FP& t, FP& dt,
              Eigen::Ref<Eigen::VectorX<FP>> ytp1) const override
    {
        assert(0 <= this->get_dt_min());
        assert(this->get_dt_min() <= this->get_dt_max());

        if (dt < this->get_dt_min() || dt > this->get_dt_max()) {
            mio::log_warning("IntegratorCore: Restricting given step size dt = {} to [{}, {}].", dt, this->get_dt_min(),
                             this->get_dt_max());
        }

        using std::max;
        using std::min;
        using std::pow;

        dt = min<FP>(dt, this->get_dt_max());

        FP t_eval; // shifted time for evaluating yt
        FP dt_new; // updated dt

        bool converged     = false; // carry for convergence criterion
        bool dt_is_invalid = false;

        if (m_yt_eval.size() != yt.size()) {
            m_yt_eval.resize(yt.size());
            m_kt_values.resize(yt.size(), m_tab_final.entries_low.size());
        }

        m_yt_eval = yt;

        while (!converged && !dt_is_invalid) {
            if (dt < this->get_dt_min()) {
                dt_is_invalid = true;
                dt            = this->get_dt_min();
            }
            // compute first column of kt, i.e. kt_0 for each y in yt_eval
            f(m_yt_eval, t, m_kt_values.col(0));

            for (Eigen::Index i = 1; i < m_kt_values.cols(); i++) {
                // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
                t_eval = t;
                t_eval += m_tab.entries[i - 1][0] *
                          dt; // t_eval = t + c_i * h // note: line zero of Butcher tableau not stored in array
                // use ytp1 as temporary storage for evaluating m_kt_values[i]
                ytp1 = m_yt_eval;
                for (Eigen::Index k = 1; k < m_tab.entries[i - 1].size(); k++) {
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

            if (converged || dt_is_invalid) {
                // if sufficiently exact, return ytp1, which currently contains the lower order approximation
                // (higher order is not always higher accuracy)
                t += dt; // this is the t where ytp1 belongs to
            }
            // else: repeat the calculation above (with updated dt)

            // compute new value for dt
            // converged implies eps/error_estimate >= 1, so dt will be increased for the next step
            // hence !converged implies 0 < eps/error_estimate < 1, strictly decreasing dt
            dt_new = dt * pow((m_eps / m_error_estimate).minCoeff(), (1. / (m_tab_final.entries_low.size() - 1)));
            // safety factor for more conservative step increases,
            // and to avoid dt_new -> dt for step decreases when |error_estimate - eps| -> 0
            dt_new *= 0.9;
            // check if updated dt stays within desired bounds and update dt for next step
            dt = min<FP>(dt_new, this->get_dt_max());
        }
        dt = max<FP>(dt, this->get_dt_min());
        // return 'converged' in favor of '!dt_is_invalid', as these values only differ if step sizing failed,
        // but the step with size dt_min was accepted.
        return converged;
    }

protected:
    Tableau<FP> m_tab;
    TableauFinal<FP> m_tab_final;
    FP m_abs_tol, m_rel_tol;
    mutable Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> m_kt_values;
    mutable Eigen::VectorX<FP> m_yt_eval;

private:
    mutable Eigen::Array<FP, Eigen::Dynamic, Eigen::Dynamic> m_eps,
        m_error_estimate; // tolerance and estimate used for time step adaption
};

} // namespace mio

#endif // ADAPT_RK_H_
