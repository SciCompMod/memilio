#ifndef EULER_H
#define EULER_H

#include <vector>
#include <epidemiology/integrator.h>

/**
 * @brief Simple explicit euler integration y(t+1) = y(t) + h*f(t,y) for ODE y'(t) = f(t,y)
 */
template <typename T> class EulerIntegrator
{
public:
    /**
     * @brief Setting up the integrator
     * @param func The right hand side of the ODE
     */
    EulerIntegrator(DerivFunction<T> func)
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
    bool step(std::vector<T> const& yt, T& t, T& dt, std::vector<T>& ytp1) const
    {
        // we are misusing the next step y as temporary space to store the derivative
        f(yt, t, ytp1);
        for (size_t i = 0; i < yt.size(); i++) {
            ytp1[i] = yt[i] + dt * ytp1[i];
        }
        t += dt;
        return true;
    }

private:
    DerivFunction<T> f;
};

#endif // EULER_H
