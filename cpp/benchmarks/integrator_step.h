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
#ifndef INTEGRATOR_STEP_H_
#define INTEGRATOR_STEP_H_

#include "memilio/io/json_serializer.h"
#include "memilio/utils/logging.h"

#include "benchmark/benchmark.h"

namespace mio
{
namespace benchmark
{
/// @brief parameters for integrator-step benchmark
struct IntegratorStepConfig {
    int num_agegroups;
    double t_init, dt_init, abs_tol, rel_tol, dt_min, dt_max;
    Eigen::VectorXd yt, ytp1;
    /**
         * @brief creates configuration with default parameters for a secir model
         * @return configuration for integrator-step benchmark
         */
    static IntegratorStepConfig initialize()
    {
        const double vals[8] = {6377.873644, 35.249156, 30.029611,   182.145865,
                                66.153059,   79.530621, 3069.383604, 159.634440};
        return IntegratorStepConfig{1,
                                    50,
                                    1,
                                    1e-10,
                                    1e-5,
                                    std::numeric_limits<double>::min(),
                                    std::numeric_limits<double>::max(),
                                    Eigen::Matrix<double, 8, 1>(vals),
                                    Eigen::VectorXd::Zero(8)};
    }
    /**
         * @brief reads configuration from json file
         * @param path the path of the configfile
         * @return configuration for integrator-step benchmark
         */
    static IntegratorStepConfig initialize(std::string path)
    {
        auto result = mio::read_json(path, mio::Tag<IntegratorStepConfig>{});
        if (!result) { // failed to read config
            mio::log(mio::LogLevel::critical, result.error().formatted_message());
            abort();
        }
        return result.value();
    }
    /// @brief function implementing mio::deserialize, used by read_json in initialize
    template <class IOContext>
    static mio::IOResult<IntegratorStepConfig> deserialize(IOContext& io)
    {
        auto obj              = io.expect_object("integrator_step");
        auto num_agegroups_io = obj.expect_element("num_agegroups", mio::Tag<int>{});
        auto t_init_io        = obj.expect_element("t_init", mio::Tag<double>{});
        auto dt_init_io       = obj.expect_element("dt_init", mio::Tag<double>{});
        auto abs_tol_io       = obj.expect_element("abs_tol", mio::Tag<double>{});
        auto rel_tol_io       = obj.expect_element("rel_tol", mio::Tag<double>{});
        auto dt_min_io        = obj.expect_element("dt_min", mio::Tag<double>{});
        auto dt_max_io        = obj.expect_element("dt_max", mio::Tag<double>{});
        auto yt_io            = obj.expect_list("yt", mio::Tag<double>{});
        return mio::apply(
            io,
            [](auto&& num_agegroups, auto&& t_init, auto&& dt_init, auto&& abs_tol, auto&& rel_tol, auto&& dt_min,
               auto&& dt_max, auto&& yt) {
                IntegratorStepConfig cfg{num_agegroups,
                                         t_init,
                                         dt_init,
                                         abs_tol,
                                         rel_tol,
                                         dt_min,
                                         dt_max,
                                         Eigen::VectorXd::Zero(yt.size()),
                                         Eigen::VectorXd::Zero(yt.size())};
                for (size_t i = 0; i < yt.size(); i++) {
                    cfg.yt[i] = yt[i];
                }
                return cfg;
            },
            num_agegroups_io, t_init_io, dt_init_io, abs_tol_io, rel_tol_io, dt_min_io, dt_max_io, yt_io);
    }
};
} // namespace benchmark

} // namespace mio

#endif
