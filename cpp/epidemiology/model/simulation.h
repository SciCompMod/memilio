#ifndef SIMULATION_H
#define SIMULATION_H

#include "epidemiology/utils/ScalarType.h"
#include "epidemiology/utils/time_series.h"
#include "epidemiology/math/adapt_rk.h"
#include "epidemiology/math/euler.h"

namespace epi
{

/**
 * @brief A class for the simulation of a compartment model.
 */
template <class Model>
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
        : m_integratorCore(std::make_shared<RKIntegratorCore>())
        , m_integrator([model](auto&& y, auto&& t, auto&& dydt) { model.eval_right_hand_side(y, t, dydt); }, t0,
                       model.get_initial_values(), dt, m_integratorCore)
        , m_model(model)
    {}

    /**
     * @brief set the core integrator used in the simulation
     */
    void set_integrator(std::shared_ptr<IntegratorCore> integrator)
    {
        m_integratorCore = std::move(integrator);
        m_integrator.set_integrator(m_integratorCore);
    }

    /**
     * @brief get_integrator
     * @return pointer to the core integrator used in the simulation
     */
    IntegratorCore* get_integrator()
    {
        return m_integratorCore.get();
    }

    /**
     * @brief get_integrator
     * @return pointer to the core integrator used in the simulation
     */
    IntegratorCore const* get_integrator() const
    {
        return m_integratorCore.get();
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
template <class Model>
TimeSeries<ScalarType> simulate(double t0, double tmax, double dt, Model const& model,
                                std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    model.check_constraints();
    Simulation<Model> sim(model, t0, dt);
    if (integrator) {
        sim.set_integrator(integrator);
    }
    sim.advance(tmax);
    return sim.get_result();
}

} // namespace epi

#endif // POPULATIONS_H
