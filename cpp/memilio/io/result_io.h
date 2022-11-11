/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kuehn
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
#ifndef MEMILIO_IO_RESULT_IO_H
#define MEMILIO_IO_RESULT_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_HDF5

#include "memilio/math/eigen_util.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/io.h"

namespace mio
{

/**
 * @brief save results of a graph simulation to h5 file
 * @param result simulation results per node of the graph.
 * @param ids identifier of each node of the graph. 
 * @param num_groups number of groups in the results.
 * @param filename name of file
 */
IOResult<void> save_result(const std::vector<TimeSeries<double>>& result, const std::vector<int>& ids, int num_groups,
                           const std::string& filename);

class SimulationResult
{
public:
    /**
     * @brief Standard constructor of SimulationResult
     * @param num_groups Number of groups or subpopulations in the simulation.
     * @param num_infectionstates Number of infection states in the considered simulation.
     */
    SimulationResult(int num_groups, int num_infectionstates)
        : m_groups(num_groups * num_infectionstates)
        , m_totals(num_infectionstates)
    {
    }

    /**
     * @brief Constructor of SimulationResult storing time, groups, and total sums of all groups
     * @param groups Simulation results of individual groups.
     * @param total Simulation results as the sum over all groups.
     */
    SimulationResult(const TimeSeries<double>& groups, const TimeSeries<double>& totals)
        : m_groups(groups)
        , m_totals(totals)
    {
    }

    /**
     * @brief Simulation results of individual groups.
     */
    const TimeSeries<double>& get_groups() const
    {
        return m_groups;
    }

    /**
     * @brief Simulation results of the sum over all groups.
     */
    const TimeSeries<double>& get_totals() const
    {
        return m_totals;
    }

private:
    TimeSeries<double> m_groups;
    TimeSeries<double> m_totals;
};

/**
 * @brief Read simulation result from h5 file.
 * @param filename name of the file to be read.
 */
IOResult<std::vector<SimulationResult>> read_result(const std::string& filename);

} // namespace mio

#endif // MEMILIO_HAS_HDF5

#endif // MEMILIO_IO_RESULT_IO_H
