/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Martin J Kuehn, Anna Wendler, Lena Ploetzke
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
#ifndef IDE_SECIR_SIMULATION_H
#define IDE_SECIR_SIMULATION_H

#include "ide_secir/parameters.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/model_ide.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include <memory>
#include <cstdio>
#include <iostream>

namespace mio
{
namespace isecir
{

/**
 * run the simulation in discrete steps and report results.
 */
class Simulation
{

public:
    /**
     * @brief setup the Simulation for an IDE model.
     * @param[in] model An instance of the IDE model.
     * @param[in] dt Step size of numerical solver.
     */
    Simulation(Model const& model, ScalarType dt = 0.1)
        : m_model(std::make_unique<Model>(model))
        , m_dt(dt)
    {
    }

    /** 
     * Run the simulation from the current time to tmax.
     * @param tmax Time to stop.
     */
    void advance(ScalarType tmax);

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all #InfectionState%s.
     * @return The result of the simulation.
     */
    TimeSeries<ScalarType>& get_result()
    {
        return m_model->m_populations;
    }

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all #InfectionState%s.
     * @return The result of the simulation.
     */
    const TimeSeries<ScalarType>& get_result() const
    {
        return m_model->m_populations;
    }

    /**
     * @brief Get the transitions between the different #InfectionState%s.
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
     * @brief get the time step of the simulation.
     * 
     */
    ScalarType get_dt()
    {
        return m_dt;
    }

private:
    std::unique_ptr<Model> m_model; ///< Unique pointer to the Model simulated.
    ScalarType m_dt; ///< Time step used for numerical computations in simulation.
};

/**
 * @brief Run a Simulation of an IDE-SECIR model.
 *
 * @param[in] tmax End time.
 * @param[in] dt Initial step size of integration.
 * @param[in] model An instance of an IDE-SECIR model.
 * @return A TimeSeries to represent the final simulation result.
 */
TimeSeries<ScalarType> simulate(double tmax, double dt, Model const& model);

} // namespace isecir
} // namespace mio

#endif //IDE_SECIR_SIMULATION_H
