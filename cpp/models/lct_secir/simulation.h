/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef LCT_SECIR_SIMULATION_H
#define LCT_SECIR_SIMULATION_H

#include "lct_secir/parameters.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/math/stepper_wrapper.h"
#include <memory>
#include <cstdio>
#include <iostream>
#include <string>

namespace mio
{
namespace lsecir
{
/**
 * @brief A class for the simulation of a LCT model.
 */
class Simulation
{
public:
    /**
     * @brief setup the simulation with an ODE solver
     * @param[in] model An instance of a LCT model
     * @param[in] t0 start time
     * @param[in] dt initial step size of integration
     */
    Simulation(Model const& model, double t0 = 0., double dt = 0.1)
        : m_integratorCore(
              std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>())
        , m_model(std::make_unique<Model>(model))
        , m_integrator(
              [&model = *m_model](auto&& y, auto&& t, auto&& dydt) {
                  model.eval_right_hand_side(y, t, dydt);
              },
              t0, m_model->get_initial_values(), dt, m_integratorCore)
    {
    }

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
     * @return reference to the core integrator used in the simulation
     */
    IntegratorCore& get_integrator()
    {
        return *m_integratorCore;
    }

    /**
     * @brief get_integrator
     * @return reference to the core integrator used in the simulation
     */
    IntegratorCore const& get_integrator() const
    {
        return *m_integratorCore;
    }

    /**
     * @brief advance simulation to tmax
     * tmax must be greater than get_result().get_last_time_point()
     * @param tmax next stopping point of simulation
     */
    Eigen::Ref<Eigen::VectorXd> advance(double tmax)
    {
        return m_integrator.advance(tmax);
    }

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all InfectionState%s.
     * @return The result of the simulation.
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
        return *m_model;
    }

    /**
     * @brief returns the simulation model used in simulation
     */
    Model& get_model()
    {
        return *m_model;
    }

    /**
     * @brief returns the next time step chosen by integrator
    */
    double get_dt() const
    {
        return m_integrator.get_dt();
    }

private:
    std::shared_ptr<IntegratorCore> m_integratorCore;
    std::unique_ptr<Model> m_model;
    OdeIntegrator m_integrator;
}; 

//TODO: In utils oder so verschieben
void print_TimeSeries(const TimeSeries<double>& result, std::string heading)
{ 
    // print compartments after simulation
    std::cout << heading << std::endl;
    for (Eigen::Index i = 0; i < result.get_num_time_points(); ++i) {
        std::cout << result.get_time(i);
        for (Eigen::Index j = 0; j < result.get_num_elements(); ++j) {
            std::cout << " | " << std::fixed << std::setprecision(8) << result[i][j];
        }
        std::cout << "\n" << std::endl;
    }
}

/**
 * @brief simulate simulates a compartmental model
 * @param[in] t0 start time
 * @param[in] tmax end time
 * @param[in] dt initial step size of integration
 * @param[in] model: An instance of a compartmental model
 * @return a TimeSeries to represent the final simulation result
 */
TimeSeries<ScalarType> simulate(double t0, double tmax, double dt, Model const& model,
                                std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    model.check_constraints();
    Simulation sim(model, t0, dt);
    if (integrator) {
        sim.set_integrator(integrator);
    }
    sim.advance(tmax);
    print_TimeSeries(sim.get_result(), model.get_heading());
    return sim.get_result();
}
} // namespace lsecir
} // namespace mio

#endif // LCT_SECIR_SIMULATION_H
