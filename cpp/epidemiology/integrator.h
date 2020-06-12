#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <memory>
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

class IntegratorCore
{
public:
    virtual ~IntegratorCore(){};

    /**
     * @brief Step of the integration with possible adaptive with
     *
     * @param[in] f Right hand side of ODE
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time step h=dt
     * @param[in,out] dt current time step h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
    virtual bool step(const DerivFunction& f, const Eigen::VectorXd& yt, double& t, double& dt,
                      Eigen::VectorXd& ytp1) const = 0;
};

/**
 * Integrate initial value problems (IVP) of ordinary differential equations (ODE) of the form y' = f(y, t), y(t0) = y0.
 */
class OdeIntegrator
{
public:
    /**
     * @brief create an integrator for a specific IVP
     * @param f rhs of the ODE
     * @param t0 initial point of independent variable t
     * @param y0 value of y at t0
     * @param dt_init initial integration step size
     * @param core implements the solution method
     */
    template <class F, class Vector>
    OdeIntegrator(F&& f, double t0, Vector&& y0, double dt_init, std::shared_ptr<IntegratorCore> core)
        : m_f(std::forward<F>(f))
        , m_t(1, t0)
        , m_y(1, std::forward<Vector>(y0))
        , m_dt(dt_init)
        , m_core(core)
    {
    }

    /**
     * @brief advance the integrator.
     * @param tmax end point. must be greater than get_t().back()
     */
    Eigen::VectorXd& advance(double tmax);

    const std::vector<double> get_t() const
    {
        return m_t;
    }

    const std::vector<Eigen::VectorXd> get_y() const
    {
        return m_y;
    }

private:
    std::shared_ptr<IntegratorCore> m_core;
    double m_dt;
    std::vector<double> m_t;
    std::vector<Eigen::VectorXd> m_y;
    DerivFunction m_f;
};

} // namespace epi

#endif // INTEGRATOR_H
