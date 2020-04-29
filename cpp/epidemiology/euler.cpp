#include <epidemiology/euler.h>

namespace epi
{

bool EulerIntegrator::step(const std::vector<double>& yt, double& t, double& dt, std::vector<double>& ytp1) const
{
    // we are misusing the next step y as temporary space to store the derivative
    f(yt, t, ytp1);
    for (size_t i = 0; i < yt.size(); i++) {
        ytp1[i] = yt[i] + dt * ytp1[i];
    }
    t += dt;
    return true;
}

} // namespace epi
