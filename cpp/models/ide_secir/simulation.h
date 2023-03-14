/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin J Kuehn
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
#ifndef EPI_ABM_SIMULATOR_H
#define EPI_ABM_SIMULATOR_H

#include "ide_secir/parameters.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "ide_secir/simulation.h"
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
 * run the simulation in discrete steps, evolve the world and report results.
 */
class Simulation
{

public:
    /**
     * @brief setup the Simulation for an IDE model.
     * @param[in] model An instance of the IDE model.
     * @param[in] t0 Start time.
     * @param[in] dt Step size of numerical solver.
     */
    Simulation(Model const& model, ScalarType t0 = 0., ScalarType dt = 0.1)
        : m_model(std::make_unique<Model>(model))
        , m_t0(t0)
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
     * Return the number of persons in all InfectionState%s.
     * @return The result of the simulation.
     */
    TimeSeries<double>& get_result()
    {
        return m_model->m_population;
    }

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all InfectionState%s.
     * @return The result of the simulation.
     */
    const TimeSeries<double>& get_result() const
    {
        return m_model->m_population;
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
     * @brief Print the transition part of the simulation result.
     * 
     * The TimeSeries m_transitions with initial values used for the simulation and calculated transitions by the 
     * simulation are printed. 
     */
    void print_transitions() const;

    /**
     * @brief Print the simulated numbers of individuals in each compartment for each time step.
     * 
     * The TimeSeries m_population with simulated numbers of individuals in each compartment for each time step are printed. 
     */
    void print_compartments() const;

private:
    std::unique_ptr<Model> m_model; ///< Unique pointer to the Model simulated.
    ScalarType m_t0; ///< Start time used for simulation.
    ScalarType m_dt; ///< Time step used for numerical computations in simulation.
};

/**
 * @brief simulate simulates a compartmental model
 * @param[in] t0 start time
 * @param[in] tmax end time
 * @param[in] dt initial step size of integration
 * @param[in] model: An instance of a compartmental model
 * @return a TimeSeries to represent the final simulation result
 * @tparam Model a compartment model type
 * @tparam Sim a simulation type that can simulate the model.
 */
TimeSeries<ScalarType> simulate(double t0, double tmax, double dt, Model const& model);

} // namespace isecir
} // namespace mio

#endif
