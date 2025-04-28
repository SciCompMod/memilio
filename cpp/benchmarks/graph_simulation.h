/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef MIO_GRAPH_SIMULATION_CONFIG_H
#define MIO_GRAPH_SIMULATION_CONFIG_H

#include "memilio/io/json_serializer.h"
#include "memilio/utils/logging.h"

#include "benchmark/benchmark.h"

namespace mio
{
namespace benchmark
{
/// @brief parameters for simulation benchmark
struct GraphConfig {
    int num_agegroups, num_regions;
    double t0, t_max, dt;
    /**
     * @brief creates configuration with default parameters for a secir model
     * @param num_agegroups number of agegroups
     * @return configuration for simulation benchmark
     */
    static GraphConfig initialize(int num_agegroups = 10)
    {
        return GraphConfig{
            num_agegroups, 400, 0, 100, 0.5,
        };
    }
    /**
     * @brief reads configuration from json file
     * @param path the path of the configfile
     * @return configuration for graph simulation benchmark
     */
    static GraphConfig initialize(std::string path)
    {
        auto result = mio::read_json(path, mio::Tag<GraphConfig>{});
        if (!result) { // failed to read config
            mio::log(mio::LogLevel::critical, result.error().formatted_message());
            abort();
        }
        return result.value();
    }
    /// @brief function implementing mio::deserialize, used by read_json in initialize
    template <class IOContext>
    static mio::IOResult<GraphConfig> deserialize(IOContext& io)
    {
        auto obj              = io.expect_object("bench_graph_simulation");
        auto num_agegroups_io = obj.expect_element("num_agegroups", mio::Tag<int>{});
        auto num_regions_io   = obj.expect_element("num_regions", mio::Tag<int>{});
        auto t_io             = obj.expect_element("t0", mio::Tag<double>{});
        auto t_max_io         = obj.expect_element("t_max", mio::Tag<double>{});
        auto dt_io            = obj.expect_element("dt", mio::Tag<double>{});
        return mio::apply(
            io,
            [](auto&& num_agegroups, auto&& num_regions, auto&& t0, auto&& t_max, auto&& dt) {
                return GraphConfig{num_agegroups, num_regions, t0, t_max, dt};
            },
            num_agegroups_io, num_regions_io, t_io, t_max_io, dt_io);
    }
};
} // namespace benchmark

} // namespace mio

#endif // MIO_GRAPH_SIMULATION_CONFIG_H
