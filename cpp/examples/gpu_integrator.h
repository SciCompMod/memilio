#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

struct tableau{
    tableau() : entries_high(6), entries_low(6), entries(5) {
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

    std::vector<double> entries_high, entries_low;
    std::vector<std::vector<double>> entries;
};


void print(std::vector<double>& vec) {
    for (auto& val : vec) {
        printf("%e ", val);
    }
    std::cout << "\n";
}

void print(std::vector<std::vector<double>>& vec) {
    for (auto& val : vec) {
        print(val);
    }
}

struct Monstrosity {

double abs_tol, rel_tol, dt_min, dt_max;

std::vector<double> m_yt_eval, m_error_estimate;
std::vector<std::vector<double>> m_kt_values; // col major!

tableau tab;

using DerivFunction = void(*)(const std::vector<double>& y, double t, std::vector<double>& dydt);

bool step(DerivFunction f, const std::vector<double>& yt, double& t, double& dt,
          std::vector<double>& ytp1) 
{
    assert(0 <= dt_min);
    assert(dt_min <= dt_max);

    if (dt < dt_min || dt > dt_max) {
        // mio::log_warning("IntegratorCore: Restricting given step size dt = {} to [{}, {}].", dt, dt_min, dt_max);
    }

    dt = std::min(dt, dt_max);

    double t_eval; // shifted time for evaluating yt
    double dt_new; // updated dt

    bool converged     = false; // carry for convergence criterion
    bool dt_is_invalid = false;

    // if (m_yt_eval.size() != yt.size()) {
    //     m_yt_eval.resize(yt.size());
    //     m_kt_values.resize(yt.size(), tab.entries_low.size());
    // }

    for (size_t j = 1; j < yt.size(); j++) {
        m_yt_eval[j] = yt[j];
    }

    while (!converged && !dt_is_invalid) {
        if (dt < dt_min) {
            dt_is_invalid = true;
            dt            = dt_min;
        }
        // std::cout << "---- step t:" << t << " dt:" << dt << "\n";
        // std::cin.ignore();
        // compute first column of kt, i.e. kt_0 for each y in yt_eval
        f(m_yt_eval, t, m_kt_values[0]);

#pragma acc loop
        for (size_t i = 1; i < m_kt_values.size(); i++) {
            // we first compute k_n1 for each y_j, then k_n2 for each y_j, etc.
            t_eval = t;
            t_eval += tab.entries[i - 1][0] *
                      dt; // t_eval = t + c_i * h // note: line zero of Butcher tableau not stored in array
            // use ytp1 as temporary storage for evaluating m_kt_values[i]
            ytp1 = m_yt_eval;
            for (size_t k = 1; k < tab.entries[i - 1].size(); k++) {
                for (size_t j = 1; j < yt.size(); j++) {
                    ytp1[j] += (dt * tab.entries[i - 1][k]) * m_kt_values[k - 1][j];
                }
            }
            // get the derivatives, i.e., compute kt_i for all y in ytp1: kt_i = f(t_eval, ytp1_low)
            f(ytp1, t_eval, m_kt_values[i]);
        }

        // for (int i = 0; i < 6; i++) {
        //     std::cout << "eval (" << i << "): ";
        //     print(m_kt_values[i]); 
        // }
        // calculate low order estimate
#pragma acc parallel loop
        for (size_t i = 0; i < yt.size(); i++) {
            ytp1[i] = m_yt_eval[i];
            for (size_t j = 0; j < m_kt_values.size(); j++) {
                ytp1[i] += (dt * (m_kt_values[j][i] * tab.entries_low[j]));
            }
        }
        // std::cout << "low: "; print(ytp1);
        // truncation error estimate: yt_low - yt_high = O(h^(p+1)) where p = order of convergence
// #pragma acc parallel loop
        for (size_t i = 0; i < yt.size(); i++) {
            m_error_estimate[i] = 0;
            for (size_t j = 0; j < m_kt_values.size(); j++) {
                m_error_estimate[i] += dt * m_kt_values[j][i] * (tab.entries_high[j] - tab.entries_low[j]);
            }
            m_error_estimate[i] = std::abs(m_error_estimate[i]);
        }

        // std::cout << "tab diff: ";
        // for (int j = 0; j < 6; j++) {
        //     std::cout << tab.entries_high[j] - tab.entries_low[j] << " ";
        // }
        // std::cout << "\n";

        // std::cout << "kt * tab: ";
        // for (int i = 0; i < yt.size(); i++) {
        //     double x = 0;
        //     for (int j = 0; j < 6; j++) {
        //         x += m_kt_values[j][i] * (tab.entries_high[j] - tab.entries_low[j]);
        //     }
        //     std::cout << x << " ";
        // }
        // std::cout << "\n";

        // std::cout << "err: "; print(m_error_estimate);
        // calculate mixed tolerance
        
        double min_coeff = std::numeric_limits<double>::infinity();
        // std::cout << "con: ";
        // #pragma acc parallel loop
        for (size_t i = 0; i < yt.size(); i++) {
            double tmp = (abs_tol + std::abs(ytp1[i]) * rel_tol) / m_error_estimate[i] ;
            // std::cout << tmp << " ";
            min_coeff = std::min(min_coeff, tmp);
        }
        converged = (min_coeff >= 1);
        // std::cout << "\nmin con: " << min_coeff << "\n";
        // converged = (min_coeff <= 1); // convergence criterion

    
        if (converged || dt_is_invalid) {
            // if sufficiently exact, return ytp1, which currently contains the lower order approximation
            // (higher order is not always higher accuracy)
            t += dt; // this is the t where ytp1 belongs to
        }
        // else: repeat the calculation above (with updated dt)

        // compute new value for dt
        // converged implies eps/error_estimate >= 1, so dt will be increased for the next step
        // hence !converged implies 0 < eps/error_estimate < 1, strictly decreasing dt
        dt_new = dt * std::pow(min_coeff, (1. / (tab.entries_low.size() - 1)));
        // safety factor for more conservative step increases,
        // and to avoid dt_new -> dt for step decreases when |error_estimate - eps| -> 0
        dt_new *= 0.9;
        // std::cout << "dt: " << dt << " dt_new: " << dt_new << "\n";
        // check if updated dt stays within desired bounds and update dt for next step
        dt = std::min(dt_new, dt_max);
    }
    dt = std::max(dt, dt_min);
    // return 'converged' in favor of '!dt_is_invalid', as these values only differ if step sizing failed,
    // but the step with size dt_min was accepted.
    return converged;
}

};
