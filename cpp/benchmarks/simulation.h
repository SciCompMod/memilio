/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef ODE_SIMULATION_CONFIG_H
#define ODE_SIMULATION_CONFIG_H

#include "memilio/io/json_serializer.h"
#include "memilio/utils/logging.h"

#include "benchmark/benchmark.h"

namespace mio
{
namespace benchmark
{
/// @brief parameters for simulation benchmark
struct SimulationConfig {
    int num_agegroups;
    double t0, t_max, dt, abs_tol, rel_tol, dt_min, dt_max;
    /**
         * @brief creates configuration with default parameters for a secir model
         * @param num_agegroups number of agegroups
         * @return configuration for simulation benchmark
         */
    static SimulationConfig initialize(int num_agegroups = 10)
    {
        return SimulationConfig{num_agegroups,
                                0,
                                100,
                                0.5,
                                1e-10,
                                1e-5,
                                std::numeric_limits<double>::min(),
                                std::numeric_limits<double>::max()};
    }
    /**
         * @brief reads configuration from json file
         * @param path the path of the configfile
         * @return configuration for simulation benchmark
         */
    static SimulationConfig initialize(std::string path)
    {
        auto result = mio::read_json(path, mio::Tag<SimulationConfig>{});
        if (!result) { // failed to read config
            mio::log(mio::LogLevel::critical, result.error().formatted_message());
            abort();
        }
        return result.value();
    }
    /// @brief function implementing mio::deserialize, used by read_json in initialize
    template <class IOContext>
    static mio::IOResult<SimulationConfig> deserialize(IOContext& io)
    {
        auto obj              = io.expect_object("bench_simulation");
        auto num_agegroups_io = obj.expect_element("num_agegroups", mio::Tag<int>{});
        auto t_io             = obj.expect_element("t0", mio::Tag<double>{});
        auto t_max_io         = obj.expect_element("t_max", mio::Tag<double>{});
        auto dt_io            = obj.expect_element("dt", mio::Tag<double>{});
        auto abs_tol_io       = obj.expect_element("abs_tol", mio::Tag<double>{});
        auto rel_tol_io       = obj.expect_element("rel_tol", mio::Tag<double>{});
        auto dt_min_io        = obj.expect_element("dt_min", mio::Tag<double>{});
        auto dt_max_io        = obj.expect_element("dt_max", mio::Tag<double>{});
        return mio::apply(
            io,
            [](auto&& num_agegroups, auto&& t0, auto&& t_max, auto&& dt, auto&& abs_tol, auto&& rel_tol, auto&& dt_min,
               auto&& dt_max) {
                return SimulationConfig{num_agegroups, t0, t_max, dt, abs_tol, rel_tol, dt_min, dt_max};
            },
            num_agegroups_io, t_io, t_max_io, dt_io, abs_tol_io, rel_tol_io, dt_min_io, dt_max_io);
    }
};
} // namespace benchmark

} // namespace mio

#endif // ODE_SIMULATION_CONFIG_H
