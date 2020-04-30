#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <functional>

/**
 * Function template to be integrated
 */
using DerivFunction = std::function<void(std::vector<double> const& y, const double t, std::vector<double>& dydt)>;

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
    virtual bool step(std::vector<double> const& yt, double& t, double& dt, std::vector<double>& ytp1) const = 0;

protected:
    DerivFunction f;
};

#endif // INTEGRATOR_H
