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
     * @brief setup the simulation with an ODE solver
     * @param[in] model An instance of a LCT model
     * @param[in] t0 start time
     * @param[in] dt initial step size of integration
     */
    Simulation(Model const& model, ScalarType t0 = 0., ScalarType dt = 0.1);

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
    Eigen::Ref<Eigen::VectorXd> advance(ScalarType tmax)
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
    ScalarType get_dt() const
    {
        return m_integrator.get_dt();
    }

private:
    std::shared_ptr<IntegratorCore> m_integratorCore;
    std::unique_ptr<Model> m_model;
    OdeIntegrator m_integrator;
};

void print_TimeSeries(const TimeSeries<ScalarType>& result, std::string heading);

/**
 * @brief simulate simulates a compartmental model
 * @param[in] t0 start time
 * @param[in] tmax end time
 * @param[in] dt initial step size of integration
 * @param[in] model: An instance of a compartmental model
 * @return a TimeSeries to represent the final simulation result
 */
TimeSeries<ScalarType> simulate(ScalarType t0, ScalarType tmax, ScalarType dt, Model const& model,
                                std::shared_ptr<IntegratorCore> integrator = nullptr);

} // namespace lsecir
} // namespace mio

#endif // LCT_SECIR_SIMULATION_H
