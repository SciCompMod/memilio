#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <functional>
#include <algorithm>
#include <Eigen/Core>

namespace epi
{

/**
 * Function template to be integrated
 */
using DerivFunction = std::function<void(const Eigen::VectorXd& y, double t, Eigen::VectorXd& dydt)>;

class IntegratorBase
{
public:
    IntegratorBase(DerivFunction func)
        : f(func)
    {
    }

    /**
     * @brief Step of the integration with possible adaptive with
     *
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time step h=dt
     * @param[in,out] dt current time step h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
    virtual bool step(const Eigen::VectorXd& yt, double& t, double& dt, Eigen::VectorXd& ytp1) const = 0;

protected:
    DerivFunction f;
};

/**
 * @brief Integrates an ODE from t0 to tmax with initian step size dt and given Integrator
 *
 * y has to be aleady allocated with the initial solution at t0
 *
 * @param t0
 * @param tmax
 * @param dt
 * @param integrator
 * @param y
 * @return Array of t with same size as y
 */
std::vector<double> ode_integrate(double t0, double tmax, double dt, const IntegratorBase& integrator, std::vector<Eigen::VectorXd>& yt);

} // namespace epi

#endif // INTEGRATOR_H
