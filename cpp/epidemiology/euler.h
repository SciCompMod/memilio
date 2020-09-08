#ifndef EULER_H
#define EULER_H

#include <vector>
#include <epidemiology/integrator.h>

namespace epi
{

class SecirParams; // forward declaration such that secir.h is not required here

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
    bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt, Eigen::Ref<Eigen::VectorXd> ytp1) const override;
};


/**
 * @brief Implicit Euler integration (not generalized, adapted to SECIHURD-model)
 */
class ImplicitEulerIntegratorCore : public IntegratorCore
{
public:
    /**
     * @brief Setting up the implicit Euler integrator
     * @param dt_min Mininum time step
     * @param dt_max Maximum time step
     * @param params Paramters of the SECIR/SECIHURD model
     */
    ImplicitEulerIntegratorCore(double dt_min, double dt_max, SecirParams const& params);

    /**
     * @brief Fixed step width of the time implicit Euler time integration scheme
     *
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time step h=dt
     * @param[in,out] dt current time step h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
    bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt, Eigen::Ref<Eigen::VectorXd> ytp1) const override;

    SecirParams const& get_secir_params() const
    {
        return m_params;
    }

    /// @param tol the required absolute tolerance for the comparison with the Fehlberg approximation (actually not really required but used in SecirSimulation constructor)
    void set_abs_tolerance(double tol)
    {
        m_abs_tol = tol;
    }

    /// @param tol the required relative tolerance for the comparison with the Fehlberg approximation (actually not really required but used in SecirSimulation constructor)
    void set_rel_tolerance(double tol)
    {
        m_rel_tol = tol;
    }


private:
    double m_abs_tol, m_rel_tol;
    double m_dt_min, m_dt_max;
    const SecirParams &m_params;
};

} // namespace epi

#endif // EULER_H
