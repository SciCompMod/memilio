#ifndef ARK_H
#define ARK_H

#include "epidemiology/math/integrator.h"

#include <cstdio>
#include <vector>

namespace epi
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
class Tableau
{
public:
    std::vector<std::vector<double>> entries;

    /**
     * @brief default is Runge-Kutta-Fehlberg4(5) tableau
     */
    Tableau();
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
class TableauFinal
{
public:
    std::vector<double> entries_low;
    std::vector<double> entries_high;

    /**
     * @brief default is Runge-Kutta-Fehlberg4(5) tableau
     */
    TableauFinal();
};

/**
 * @brief Two scheme Runge-Kutta numerical integrator with adaptive step width
 *
 * This method integrates a system of ODEs
 */
class RKIntegratorCore : public IntegratorCore
{
public:
    /**
     * @brief Setting up the integrator
     * @param func The right hand side of the ODE
     * @param dt_min Mininum time step
     * @param dt_max Maximum time step
     */
    RKIntegratorCore()
        : m_abs_tol(1e-10)
        , m_rel_tol(1e-5)
        , m_dt_min(std::numeric_limits<double>::min())
        , m_dt_max(std::numeric_limits<double>::max())
    {
    }

    /// @param tol the required absolute tolerance for the comparison with the Fehlberg approximation
    void set_abs_tolerance(double tol)
    {
        m_abs_tol = tol;
    }

    /// @param tol the required relative tolerance for the comparison with the Fehlberg approximation
    void set_rel_tolerance(double tol)
    {
        m_rel_tol = tol;
    }

    /// sets the minimum step size (not used)
    void set_dt_min(double dt_min) override
    {
        m_dt_min = dt_min;
    }

    /// sets the minimum step size (not used)
    void set_dt_max(double dt_max) override
    {
        m_dt_max = dt_max;
    }

    // Allow setting different RK tablea schemes
    void set_tableaus(const Tableau& tab, const TableauFinal& final_tab)
    {
        m_tab       = tab;
        m_tab_final = final_tab;
    }

    /**
    * Adaptive step width of the integration
    * This method integrates a system of ODEs
    * @param[in] yt value of y at t, y(t)
    * @param[in,out] t current time step h=dt
    * @param[in,out] dt current time step h=dt
    * @param[out] ytp1 approximated value y(t+1)
    */
    bool step(const DerivFunction& f, Eigen::Ref<Eigen::VectorXd const> yt, double& t, double& dt,
              Eigen::Ref<Eigen::VectorXd> ytp1) const override;

private:
    Tableau m_tab;
    TableauFinal m_tab_final;
    double m_abs_tol, m_rel_tol;
    double m_dt_min, m_dt_max;
};

} // namespace epi

#endif // ARK_H
