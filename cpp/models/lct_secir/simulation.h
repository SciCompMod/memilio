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
#include <memory>
#include <cstdio>
#include <iostream>

namespace mio
{
namespace lsecir
{

/**
 * run the simulation in discrete steps and report results.
 */
class Simulation
{

public:
    /**
     * @brief setup the Simulation for an LCT model.
     * @param[in] model An instance of the LCT model.
     * @param[in] t0 Start time.
     * @param[in] dt Step size of numerical solver.
     */
    Simulation(Model const& model, ScalarType t0 = 0., ScalarType dt = 0.1)
        : m_model(std::make_unique<Model>(model))
        , m_t0(t0)
        , m_dt(dt),
        m_integratorCore(
              std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>()),
        m_integrator(
              [&model = *m_model](auto&& y, auto&& t, auto&& dydt) {
                  model.eval_right_hand_side(y, y, t, dydt);
              },
              t0, m_model->get_initial_values(), dt, m_integratorCore)
    {
    }

    /*
     * @brief set the core integrator used in the simulation
     */
    void set_integrator(std::shared_ptr<IntegratorCore> integrator)
    {
        m_integratorCore = std::move(integrator);
        m_integrator.set_integrator(m_integratorCore);
    }

    /** 
     * Run the simulation from the current time to tmax.
     * @param tmax Time to stop.
     */
    void advance(ScalarType tmax);

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all InfectionState%s.
     * @return The result of the simulation.
     */
    TimeSeries<double>& get_result()
    {
        return m_model->m_populations;
    }

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all InfectionState%s.
     * @return The result of the simulation.
     */
    const TimeSeries<double>& get_result() const
    {
        return m_model->m_populations;
    }

    /**
     * @brief Get the transitions between the different InfectionState%s.
     * 
     * @return TimeSeries with stored transitions calculated in the simulation.
     */
    TimeSeries<ScalarType> const& get_transitions()
    {
        return m_model->m_transitions;
    }

    /**
     * @brief returns the simulation model used in simulation.
     */
    const Model& get_model() const
    {
        return *m_model;
    }

    /**
     * @brief returns the simulation model used in simulation.
     */
    Model& get_model()
    {
        return *m_model;
    }

    /**
     * @brief get the starting time of the simulation.
     * 
     */
    ScalarType get_t0()
    {
        return m_t0;
    }

    /**
     * @brief get the time step of the simulation.
     * 
     */
    ScalarType get_dt()
    {
        return m_dt;
    }

    /**
     * @brief Print the transition part of the simulation result.
     * 
     * The TimeSeries m_transitions with initial values used for the simulation and calculated transitions by the 
     * simulation are printed. 
     */
    void print_transitions() const;

    /**
     * @brief Print the simulated numbers of individuals in each compartment for each time step.
     * 
     * The TimeSeries m_populations with simulated numbers of individuals in each compartment for each time step are printed. 
     */
    void print_compartments() const;

private:
    std::unique_ptr<Model> m_model; ///< Unique pointer to the Model simulated.
    ScalarType m_t0; ///< Start time used for simulation.
    ScalarType m_dt; ///< Time step used for numerical computations in simulation.
};

/**
 * @brief simulates a compartmental model
 * @param[in] t0 start time
 * @param[in] tmax end time
 * @param[in] dt initial step size of integration
 * @param[in] model an instance of a compartmental model
 * @return a TimeSeries to represent the final simulation result
 */
TimeSeries<ScalarType> simulate(double t0, double tmax, double dt, Model const& model);

} // namespace isecir
} // namespace mio

#endif //LCT_SECIR_SIMULATION_H
