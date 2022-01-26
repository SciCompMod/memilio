/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele
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
#ifndef SECIR_RESULT_IO_H
#define SECIR_RESULT_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_HDF5

#include "secir_vaccine/secir.h"
#include "memilio/math/eigen_util.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/io.h"

namespace mio
{
namespace vaccinated
{

    /**
 * @brief save secir simulation result to h5 file
 * @param times Vector of timesteps used during simulation
 * @param secir Results of secir simulation
 * @param filename name of file
 */
    IOResult<void> save_result(const std::vector<TimeSeries<double>>& result, const std::vector<int>& ids,
                               const std::string& filename);

    class SecirSimulationResult
    {
    public:
        /**
     * @brief Standard constructor of SecirSimulationResult
     */
        SecirSimulationResult(int num_groups, int num_compartments)
            : m_groups(num_groups * num_compartments)
            , m_totals(num_compartments)
        {
        }

        /**
     * @brief Constructor of SecirSimulationResult storing time, groups, and total sums of all groups
     * @param groups Simulation Results of individual groups
     * @param total Simulation Results of the sum over all groups
     */
        SecirSimulationResult(const TimeSeries<double>& groups, const TimeSeries<double>& totals)
            : m_groups(groups)
            , m_totals(totals)
        {
        }

        /**
     * @brief Simulation Results of individual groups.
     */
        const TimeSeries<double>& get_groups() const
        {
            return m_groups;
        }

        /**
     * @brief Simulation Results of the sum over all groups.
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
 * @brief read secir simulation result from h5 file
 * @param filename name of file
 * @param nb_groups number of groups used during simulation
 */
    IOResult<std::vector<SecirSimulationResult>> read_result(const std::string& filename, int nb_groups);

} // namespace vaccinated
} // namespace mio

#endif // MEMILIO_HAS_HDF5

#endif // SECIR_RESULT_IO_H
