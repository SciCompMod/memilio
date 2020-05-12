#include <epidemiology/integrator.h>

#include <epidemiology/logging.h>

namespace epi
{

std::vector<double> ode_integrate(double t0, double tmax, double dt, const IntegratorBase& integrator,
                                  std::vector<Eigen::VectorXd>& y)
{
    assert(tmax > t0);

    size_t ode_dim = y[0].size();
    size_t nb_steps = (int)(ceil((tmax - t0) / dt)); // estimated number of time steps (if equidistant)

    std::vector<double> vec_times(1, t0);
    vec_times.reserve(nb_steps + 1);
    y.reserve(nb_steps + 1);

    bool step_okay = true;

    double t      = t0;
    size_t i      = 0;
    while (std::abs((tmax - t) / (tmax - t0)) > 1e-10) {
        //we don't make timesteps too small as the error estimator of an adaptive integrator 
        //may not be able to handle it. this is very conservative and maybe unnecessary, 
        //but also unlikely to happen. may need to be reevaluated

        dt = std::min(dt, tmax - t);
        y.emplace_back(ode_dim);
        step_okay &= integrator.step(y[i], t, dt, y[i + 1]);
        vec_times.push_back(t);
        i++;
    }

    if (!step_okay) {
        log_warning("Adaptive step sizing failed.");
    } else if (std::abs((tmax - t) / (tmax - t0))  > 1e-15) {
        log_warning("Last time step too small. Could not reach tmax exactly.");
    } else {
        log_info("Adaptive step sizing successful to tolerances.");        
    }

    return vec_times;
}

} // namespace epi
