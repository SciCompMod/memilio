#ifndef ARK_H
#define ARK_H

#include <cstdio>
#include <vector>
#include <functional>

#include <epidemiology/seirParam.h>

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
template <typename T>
struct tableau {

    std::vector<std::vector<T>> entries;

    // default is Runge-Kutta-Fehlberg4(5) tableau
    tableau()
    {
        entries.resize(5);
        for(size_t i = 0; i < entries.size(); i++) {
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
template <typename T>
struct tableau_final {

    std::vector<T> entries_low;
    std::vector<T> entries_high;

// default is Runge-Kutta-Fehlberg4(5) tableau
    tableau_final()
    {
        entries_low.resize(6);
        entries_high.resize(6);

        entries_low[0] = 25 / 216.0;
        entries_low[1] = 0.0;
        entries_low[2] = 1408 / 2565.0;
        entries_low[3] = 2197 / 4104.0;
        entries_low[4] = -0.2;
        entries_low[5] = 0.0;
        entries_high[0] = 16 / 135.0;
        entries_high[1] = 0.0;
        entries_high[2] = 6656 / 12825.0;
        entries_high[3] = 28561 / 56430.0;
        entries_high[4] = -9 / 50.0;
        entries_high[5] = 2 / 55.0;
    }
};


/**
 * Function template to be integrated
 */
template <typename T>
using DerivFunction = std::function<void(std::vector<T> const &y, const T t, std::vector<T> &dydt)>;

/**
 * @brief
 * Two scheme Runge-Kutta numerical integrator with adaptive step width
 * This method integrates a system of ODEs
 */
template <typename T>
class RKIntegrator
{
public:
    /**
     * @brief Setting up the integrator
     * @param func The right hand side of the ODE
     * @param dt_min Mininum time step
     * @param dt_max Maximum time step
     */
    RKIntegrator(DerivFunction<T> func, T dt_min, T dt_max)
        : f(func), m_abs_tol(1e-10), m_rel_tol(1e-5), m_dt_min(dt_min), m_dt_max(dt_max)
    {}

    /// @param tol_abs the required absolute tolerance for the comparison with the Fehlberg approximation
    void set_abs_tolerance(T tol)
    {
        m_abs_tol = tol;
    }

    /// @param tol_rel the required relative tolerance for the comparison with the Fehlberg approximation
    void set_rel_tolerance(T tol)
    {
        m_rel_tol = tol;
    }

    // Allow setting different RK tablea schemes
    void set_tableaus(const tableau<T>& tab, const tableau_final<T>& final_tab)
    {
        m_tab = tab;
        m_tab_final = final_tab;
    }

    /**
    * Adaptive step width of the integration
    * This method integrates a system of ODEs
    * @param[in] yt value of y at t, y(t)
    * @param[in] f right hand side of ODE f(t,y)
    * @param[in] dtmin the minimum step size
    * @param[in] dtmax the maximum step size
    * @param[inout] t current time step h=dt
    * @param[inout] dt current time step h=dt
    * @param[out] ytp1 approximated value y(t+1)
    */
    bool step(std::vector<T> const& yt, T& t, T& dt, std::vector<T>& ytp1) const
    {

        T max_err = 1e10;
        T conv_crit = 1e9;

        std::vector<std::vector<T>> kt_values;

        bool failed_step_size_adapt = false;

        dt = 2 * dt;

        while(max_err > conv_crit && !failed_step_size_adapt) {
            dt = 0.5 * dt;

            kt_values.resize(0); // remove data from previous loop
            kt_values.resize(m_tab_final.entries_low.size()); // these are the k_ni per y_t, used to compute y_t+1

            for(size_t i = 0; i < kt_values.size(); i++) { // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
                kt_values[i].resize(yt.size()); // note: yt contains more than one variable since we solve a system of ODEs

                if(i == 0) { // for i==0, we have kt_i=f(t,y(t))
                    f(yt, t, kt_values[i]);
                    // printf("\n\n t = %.4f\t kn1 = %.4f", t, kt_values[i][0]);
                }
                else {

                    double t_eval = t;

                    t_eval += m_tab.entries[i - 1][0] * dt; // t_eval = t + c_i * h // note: line zero of Butcher tableau not stored in array !

                    // printf("\n t = %.4f", t_eval);

                    std::vector<T> yt_eval(yt);

                    // go through the variables of the system: S, E, I, ....
                    for(size_t j = 0; j < yt_eval.size(); j++) {
                        // go through the different kt_1, kt_2, ..., kt_i-1 and add them onto yt: y_eval = yt + h * \sum_{j=1}^{i-1} a_{i,j} kt_j
                        for(size_t k = 1; k < m_tab.entries[i - 1].size() - 1; k++) {
                            yt_eval[j] += ( dt * m_tab.entries[i - 1][k] * kt_values[k - 1][j] ); // note the shift in k and k-1 since the first column of 'tab' corresponds to 'b_i' and 'a_ij' starts with the second column
                        }

                    }

                    // get the derivatives, i.e., compute kt_i for all y at yt_eval: kt_i = f(t_eval, yt_eval)
                    f(yt_eval, t_eval, kt_values[i]);

                    // printf("\t kn%d = %.4f", i+1, kt_values[i][0]);

                }

            }

            // copy actual yt to compare both new iterates
            std::vector<T> ytp1_low(yt); // lower order approximation (e.g., order 4)
            std::vector<T> ytp1_high(yt); // higher order approximation (e.g., order 5)

            std::vector<T> err(yt.size(), 0);
            double max_val = 0;
            max_err = 0;

            for (size_t i = 0; i < yt.size(); i++) {

                for(size_t j = 0; j < kt_values.size(); j++) { // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
                    ytp1_low[i] += (dt * m_tab_final.entries_low[j] * kt_values[j][i]);
                    ytp1_high[i] += (dt * m_tab_final.entries_high[j] * kt_values[j][i]);
                }


                err[i] = 1 / dt * std::abs(ytp1_low[i] - ytp1_high[i]); // divide by h=dt since the local error is one order higher than the global one

                // printf("\n i: %lu,\t err_abs %e, err_rel %e\n", i,  err[i], max_val);
                if(err[i] > max_err) {
                    max_err = err[i];
                }
                if(max_val < ytp1_low[i]) {
                    max_val = ytp1_low[i];
                }

            }

            // printf("\n low %.8f\t high %.8f\t max_err: %.8f,\t val: %.8f,\t crit: %.8f,\t %d t: %.8f\t dt: %.8f ", ytp1_low[0], ytp1_high[0], max_err, max_val, m_abs_tol + max_val*tol_rel, max_err <= (m_abs_tol + max_val*m_rel_tol), t, dt);
            // sleep(1);
            conv_crit = m_abs_tol + max_val * m_rel_tol;

            if(max_err <= conv_crit || dt < 2 * m_dt_min + 1e-6) {
                // if sufficiently exact, take 4th order approximation (do not take 5th order : Higher order is not always higher accuracy!)
                for (size_t i = 0; i < yt.size(); i++) {
                    ytp1[i] = ytp1_low[i];
                }

                if(dt < 2 * m_dt_min + 1e-6) {
                    failed_step_size_adapt = true;
                }

                t += dt; // this is the t where ytp1 belongs to

                if(max_err <= 0.03 * conv_crit && 2 * dt < m_dt_max) { // error of doubled step size is about 32 times as large
                    dt = 2 * dt;
                }

            }

        }

        return failed_step_size_adapt;

    }

private:
    DerivFunction<T> f;
    tableau<T> m_tab;
    tableau_final<T> m_tab_final;
    T m_abs_tol, m_rel_tol;
    T m_dt_min, m_dt_max;
};

#endif // ARK_H
