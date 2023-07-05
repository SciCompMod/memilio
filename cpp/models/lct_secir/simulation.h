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

#include "lct_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/metaprogramming.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/math/eigen.h"
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
     * @brief Set up the simulation with an ODE solver.
     *
     * Per default Runge Kutta Cash Karp 54 solver is used. 
     * Use function set_integrator() to set a different integrator.
     * @param[in] model An instance of a LCT model.
     * @param[in] t0 The Start time, usually 0.
     * @param[in] dt Initial step size of integration.
     */
    Simulation(Model const& model, ScalarType t0 = 0., ScalarType dt = 0.1);

    /**
     * @brief Set the core integrator used in the simulation.
     *
     * @param[in] integrator Core integrator that should be used for simulation.
     */
    void set_integrator(std::shared_ptr<IntegratorCore> integrator)
    {
        m_integratorCore = std::move(integrator);
        m_integrator.set_integrator(m_integratorCore);
    }

    /**
     * @brief Get a reference to the used integrator.
     *
     * @return Reference to the core integrator used in the simulation
     */
    IntegratorCore& get_integrator()
    {
        return *m_integratorCore;
    }

    /**
     * @brief Get a reference to the used integrator.
     *
     * @return Reference to the core integrator used in the simulation
     */
    IntegratorCore const& get_integrator() const
    {
        return *m_integratorCore;
    }

    /**
     * @brief Advance simulation to tmax.
     *
     * tmax must be greater than get_result().get_last_time_point().
     * @param tmax Stopping point of the simulation.
     */
    Eigen::Ref<Eigen::VectorXd> advance(ScalarType tmax)
    {
        return m_integrator.advance(tmax);
    }

    /**
     * @brief Get the result of the simulation.
     *
     * Return the number of persons in all InfectionState%s (inclusive subcompartments).
     * @return The result of the simulation in form of a TimeSeries.
     */
    TimeSeries<ScalarType>& get_result()
    {
        return m_integrator.get_result();
    }

    /**
     * @brief Get the result of the simulation.
     *
     * Return the number of persons in all InfectionState%s (inclusive subcompartments).
     * @return The result of the simulation in form of a TimeSeries.
     */
    const TimeSeries<ScalarType>& get_result() const
    {
        return m_integrator.get_result();
    }

    /**
     * @brief Returns a reference to the simulation model used in simulation.
     */
    const Model& get_model() const
    {
        return *m_model;
    }

    /**
     * @brief Returns a reference to the simulation model used in simulation.
     */
    Model& get_model()
    {
        return *m_model;
    }

    /**
     * @brief Returns the next time step chosen by integrator.
     */
    ScalarType get_dt() const
    {
        return m_integrator.get_dt();
    }

private:
    std::shared_ptr<IntegratorCore> m_integratorCore; ///< InteratorCore used for Simulation.
    std::unique_ptr<Model> m_model; ///< LCT-model the simulation should be performed with.
    OdeIntegrator m_integrator; ///< OdeIntegrator used to perform simulation.
};

/**
 * @brief Prints values in a TimeSeries.
 *
 * @param[in] result A TimeSeries that should be printed.
 * @param[in] heading A heading that should be printed before the TimeSeries.
 */
void print_TimeSeries(const TimeSeries<ScalarType>& result, std::string heading = " ");

/**
 * @brief Performs a simulation with specified parameters for given model.
 *
 * @param[in] t0 Start time of the simulation.
 * @param[in] tmax End time of the simulation.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of a LCT model for which the simulation should be performed.
 * @param[in] integrator An integrator that should be used to discretize the model equations. 
 *      If default value is used the simulation will be performed with the runge_kutta_cash_karp54 method.
 * @return A TimeSeries with the result of the simulation.
 */
TimeSeries<ScalarType> simulate(ScalarType t0, ScalarType tmax, ScalarType dt, Model const& model,
                                std::shared_ptr<IntegratorCore> integrator = nullptr);

} // namespace lsecir
} // namespace mio

#endif // LCT_SECIR_SIMULATION_H
