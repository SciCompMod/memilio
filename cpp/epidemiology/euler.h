#ifndef EULER_H
#define EULER_H

#include <vector>
#include <epidemiology/integrator.h>

namespace epi
{

/**
 * @brief Simple explicit euler integration y(t+1) = y(t) + h*f(t,y) for ODE y'(t) = f(t,y)
 */
class EulerIntegrator
{
public:
    /**
     * @brief Setting up the integrator
     * @param func The right hand side of the ODE
     */
    EulerIntegrator(DerivFunction func)
        : f(func)
    {
    }

    /**
    * Adaptive step width of the integration
    * This method integrates a system of ODEs
    * @param[in] yt value of y at t, y(t)
    * @param[inout] t current time step h=dt
    * @param[inout] dt current time step h=dt
    * @param[out] ytp1 approximated value y(t+1)
    */
    bool step(std::vector<double> const& yt, double& t, double& dt, std::vector<double>& ytp1) const;

private:
    DerivFunction f;
};

} // namespace epi

#endif // EULER_H
