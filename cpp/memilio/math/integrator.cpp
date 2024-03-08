/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele
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
#include "memilio/math/integrator.h"
#include "memilio/utils/logging.h"

namespace mio
{

Eigen::Ref<Eigen::VectorXd> OdeIntegrator::advance(const DerivFunction& f, const double tmax, double& dt,
                                                   TimeSeries<double>& results)
{
    const double t0 = results.get_last_time();
    assert(tmax > t0);
    assert(dt > 0);

    const size_t nb_steps = (int)(ceil((tmax - t0) / dt)); // estimated number of time steps (if equidistant)

    results.reserve(results.get_num_time_points() + nb_steps);

    bool step_okay = true;

    double t = t0;
    size_t i = results.get_num_time_points() - 1;
    while (std::abs((tmax - t) / (tmax - t0)) > 1e-10) {
        //we don't make timesteps too small as the error estimator of an adaptive integrator
        //may not be able to handle it. this is very conservative and maybe unnecessary,
        //but also unlikely to happen. may need to be reevaluated

        auto dt_eff = std::min(dt, tmax - t);
        results.add_time_point();
        step_okay &= m_core->step(f, results[i], t, dt_eff, results[i + 1]);
        results.get_last_time() = t;

        ++i;

        if (std::abs((tmax - t) / (tmax - t0)) > 1e-10 || dt_eff > dt) {
            //store dt only if it's not the last step as it is probably smaller than required for tolerances
            //except if the step function returns a bigger step size so as to not lose efficiency
            dt = dt_eff;
        }
    }

    if (!step_okay) {
        log_warning("Adaptive step sizing failed. Forcing an integration step of size dt_min.");
    }
    else if (std::abs((tmax - t) / (tmax - t0)) > 1e-14) {
        log_warning("Last time step too small. Could not reach tmax exactly.");
    }
    else {
        log_info("Adaptive step sizing successful to tolerances.");
    }

    return results.get_last_value();
}

} // namespace mio
