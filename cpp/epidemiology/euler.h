#ifndef EULER_H
#define EULER_H

#include <vector>
#include <epidemiology/integrator.h>

namespace epi
{

/**
 * @brief Simple explicit euler integration y(t+1) = y(t) + h*f(t,y) for ODE y'(t) = f(t,y)
 */
class EulerIntegratorCore : public IntegratorCore
{
public:
    /**
     * @brief Fixed step width of the integration
     *
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time step h=dt
     * @param[in,out] dt current time step h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
    bool step(const DerivFunction& f, Eigen::VectorXd const& yt, double& t, double& dt, Eigen::VectorXd& ytp1) const override;
};

} // namespace epi

#endif // EULER_H
