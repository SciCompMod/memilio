#include <epidemiology/euler.h>

namespace epi
{

bool EulerIntegratorCore::step(const DerivFunction& f, const Eigen::VectorXd& yt, double& t, double& dt, Eigen::VectorXd& ytp1) const
{
    // we are misusing the next step y as temporary space to store the derivative
    f(yt, t, ytp1);
    ytp1 = yt + dt * ytp1;
    t += dt;
    return true;
}

} // namespace epi
