#ifndef EULER_H
#define EULER_H

#include "epidemiology/math/integrator.h"

#include <vector>

namespace epi
{

class SecirModel;

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
    bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
              Eigen::Ref<Eigen::VectorXd> ytp1) const override;
};

/**
 * @brief Implicit Euler integration (not generalized, adapted to SECIHURD-model)
 */
class ImplicitEulerIntegratorCore : public IntegratorCore
{
public:
    /**
     * @brief Setting up the implicit Euler integrator
     * @param params Paramters of the SECIR/SECIHURD model
     */
    ImplicitEulerIntegratorCore(SecirModel const& params);

    /**
     * @brief Fixed step width of the time implicit Euler time integration scheme
     *
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time step h=dt
     * @param[in,out] dt current time step h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
    bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
              Eigen::Ref<Eigen::VectorXd> ytp1) const override;

    SecirModel const& get_secir_params() const
    {
        return m_model;
    }

    /// @param tol the required absolute tolerance for the comparison with the Fehlberg approximation (actually not really required but used in SecirSimulation constructor)
    void set_abs_tolerance(double tol)
    {
        m_abs_tol = tol;
    }

private:
    const SecirModel& m_model;
    double m_abs_tol = 1e-4;
};

} // namespace epi

#endif // EULER_H
