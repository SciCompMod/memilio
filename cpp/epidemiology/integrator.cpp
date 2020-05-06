#include <epidemiology/integrator.h>

#include <epidemiology/logging.h>

namespace epi
{

std::vector<double> ode_integrate(double t0, double tmax, double dt, const IntegratorBase& integrator,
                                  std::vector<Eigen::VectorXd>& y)
{

    size_t n_params = y[0].size();
    size_t nb_steps = (int)(ceil((tmax - t0) / dt)); // estimated number of time steps (if equidistant)

    y.resize(nb_steps + 1, Eigen::VectorXd::Constant(n_params, 0));
    std::vector<double> vec_times(nb_steps + 1, 0.);

    //@TODO: get from integrator base
    double dtmin = 1e-3;

    bool step_okay = true;

    double t      = t0;
    double t_prev = t0;
    size_t i      = 0;
    vec_times[0]  = t0;
    while (t_prev < tmax) {
        if (t > tmax) { // possible for adaptive step size
            dt = tmax - t_prev;
            if (dt < 0.1 * dtmin) {
                break;
            }
        }
        t_prev = t;
        t      = std::min(t, tmax); // possible for adaptive step size

        if (i + 1 >= y.size()) {
            std::vector<Eigen::VectorXd> vecAppend(20, Eigen::VectorXd::Constant(n_params, 0));
            y.insert(y.end(), vecAppend.begin(), vecAppend.end());
            vec_times.resize(vec_times.size() + 20, 0.);
        }
        step_okay        = integrator.step(y[i], t, dt, y[i + 1]);
        vec_times[i + 1] = t;

        i++;
    }

    // cut empty elements (makes more sense for adaptive time step size)
    if (y.size() > i) {
        y.resize(i);
        vec_times.resize(i);
    }

    if (step_okay) {
        epi::log_info("Adaptive step sizing successful to tolerances.");
    }
    else {
        epi::log_info("Adaptive step sizing failed.");
    }

    return vec_times;
}

} // namespace epi
