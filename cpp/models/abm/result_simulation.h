/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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

#include "abm/common_abm_loggers.h"
#include "abm/simulation.h"

namespace mio
{
namespace abm
{

/// @brief Simulation holding its own History to provide a get_result member. Can be used for a ParameterStudy.
template <class M = Model>
class ResultSimulation : public Simulation<M>
{
public:
    using Model = M;

    /// @brief Create a simulation, moving the model.
    ResultSimulation(Model&& m, TimePoint t)
        : Simulation<Model>(t, std::move(m))
    {
    }

    /**
     * @brief Run the simulation until the given time point.
     * @param tmax Final time point for the simualtion. 
     */
    void advance(TimePoint tmax)
    {
        Simulation<Model>::advance(tmax, history);
    }

    /**
     * @brief Return the simulation result aggregated by infection states.
     */
    const mio::TimeSeries<double>& get_result() const
    {
        return get<0>(history.get_log());
    }

    mio::History<TimeSeriesWriter, LogInfectionState> history{
        Eigen::Index(InfectionState::Count)}; ///< History used to create the result TimeSeries.
};

} // namespace abm
} // namespace mio
