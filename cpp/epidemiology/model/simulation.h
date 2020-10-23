#ifndef SIMULATION_H
#define SIMULATION_H

#include "epidemiology/utils/ScalarType.h"
#include "epidemiology/utils/time_series.h"
#include "epidemiology/math/adapt_rk.h"
#include "epidemiology/math/euler.h"

namespace epi
{

using DefaultIntegratorCore = RKIntegratorCore;

/**
 * @brief A class for the simulation of a compartment model.
 */
template <class Model, class IntegratorCore = DefaultIntegratorCore>
class Simulation
{
public:
    /**
     * @brief setup the simulation with an ODE solver
     * @param[in] model: An instance of a compartmental model
     * @param[in] t0 start time
     * @param[in] dt initial step size of integration
     */
    Simulation(Model const& model, double t0 = 0., double dt = 0.1)
        : m_integratorCore(std::make_shared<IntegratorCore>())
        , m_integrator([model](auto&& y, auto&& t, auto&& dydt) { model.eval_right_hand_side(y, t, dydt); }, t0,
                       model.get_initial_values(), dt, m_integratorCore)
        , m_model(model)
    {
        m_integratorCore->set_dt_min(0.3);
        m_integratorCore->set_dt_max(1.0);
        m_integratorCore->set_rel_tolerance(1e-4);
        m_integratorCore->set_abs_tolerance(1e-1);
    }

    /**
     * @brief advance simulation to tmax
     * must be greater than get_t().back()
     * @param tmax next stopping point of simulation
     */
    Eigen::Ref<Eigen::VectorXd> advance(double tmax)
    {
        return m_integrator.advance(tmax);
    }

    /**
     * @brief get_result returns the final simulation result
     * @return a TimeSeries to represent the final simulation result
     */
    TimeSeries<ScalarType>& get_result()
    {
        return m_integrator.get_result();
    }

    /**
     * @brief get_result returns the final simulation result
     * @return a TimeSeries to represent the final simulation result
     */
    const TimeSeries<ScalarType>& get_result() const
    {
        return m_integrator.get_result();
    }

    /**
     * @brief returns the simulation model used in simulation
     */
    const Model& get_model() const
    {
        return m_model;
    }

private:
    std::shared_ptr<IntegratorCore> m_integratorCore;
    OdeIntegrator m_integrator;
    Model m_model;
}; // namespace epi

/**
 * @brief simulate simulates a compartmental model
 * @param[in] t0 start time
 * @param[in] tmax end time
 * @param[in] dt initial step size of integration
 * @param[in] model: An instance of a compartmental model
 * @return a TimeSeries to represent the final simulation result
 *
 */
template <class Model, class IntegratorCore = DefaultIntegratorCore>
TimeSeries<ScalarType> simulate(double t0, double tmax, double dt, Model const& model)
{
    model.check_constraints();
    Simulation<Model, IntegratorCore> sim(model, t0, dt);
    sim.advance(tmax);
    return sim.get_result();
}

} // namespace epi

#endif // POPULATIONS_H
