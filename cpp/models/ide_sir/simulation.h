/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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

/* This file contains the SImulation class for the model ModelMessinaExtendedDetailedInit which will be further investiagted and developed. */

#ifndef IDE_SIR_SIMULATION_H
#define IDE_SIR_SIMULATION_H

#include "ide_sir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include <cstdio>

namespace mio
{
namespace isir
{

class SimulationMessinaExtendedDetailedInit
{

public:
    /**
     * @brief setup the Simulation for an IDE model.
     * @param[in] model An instance of the IDE model.
     * @param[in] dt Step size of numerical solver.
     */
    SimulationMessinaExtendedDetailedInit(ModelMessinaExtendedDetailedInit const& model, ScalarType dt)
        : m_model(std::make_unique<ModelMessinaExtendedDetailedInit>(model))
        , m_dt(dt)
    {
        assert(m_dt > 0);
    }

    /** 
     * Run the simulation from the current time to tmax.
     * @param tmax Time to stop.
     */
    void advance(ScalarType tmax, bool backwards_fd = true, bool use_complement = false);

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all #InfectionState%s.
     * @return The result of the simulation.
     */
    TimeSeries<ScalarType> get_result()
    {
        return m_model->populations;
    }

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all #InfectionState%s.
     * @return The result of the simulation.
     */
    const TimeSeries<ScalarType>& get_result() const
    {
        return m_model->populations;
    }

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all #InfectionState%s.
     * @return The result of the simulation.
     */
    TimeSeries<ScalarType> get_flows()
    {
        return m_model->flows;
    }

    /**
     * @brief Get the result of the simulation.
     * Return the number of persons in all #InfectionState%s.
     * @return The result of the simulation.
     */
    const TimeSeries<ScalarType>& get_flows() const
    {
        return m_model->flows;
    }

    /**
     * @brief returns the simulation model used in simulation.
     */
    const ModelMessinaExtendedDetailedInit& get_model() const
    {
        return *m_model;
    }

    /**
     * @brief returns the simulation model used in simulation.
     */
    ModelMessinaExtendedDetailedInit& get_model()
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
    std::unique_ptr<ModelMessinaExtendedDetailedInit> m_model; ///< Unique pointer to the Model simulated.
    ScalarType m_dt; ///< Time step used for numerical computations in simulation.
    size_t m_max_number_iterations =
        0; ///< Get maximal number of iterations that was necessary throughout the simulation.
};

} // namespace isir
} // namespace mio

#endif //IDE_SIR_SIMULATION_H
