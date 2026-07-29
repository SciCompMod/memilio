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

#include "memilio/io/cli.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include "memilio/utils/logging.h"
#include <cmath>
#include <functional>
#include <filesystem>
#include <numbers>
#include <vector>

using namespace mio;

namespace params
{
ScalarType xleft  = 0.;
ScalarType xright = 100.;

ScalarType yleft  = 1.;
ScalarType yright = 0.;

size_t num_steps = (xright - xleft) * 1000000 + 1;

} // namespace params

ScalarType smoothercos(ScalarType normalized_time, ScalarType yleft, ScalarType yright)
{
    return 0.5 * (yleft - yright) * std::cos(std::numbers::pi_v<ScalarType> * normalized_time) + 0.5 * (yleft + yright);
}

ScalarType smoothstep_c1(ScalarType normalized_time, ScalarType yleft, ScalarType yright)
{
    return yleft + (yright - yleft) * (3. * std::pow(normalized_time, 2.) - 2. * std::pow(normalized_time, 3.));
}

ScalarType smoothstep_c2(ScalarType normalized_time, ScalarType yleft, ScalarType yright)
{
    return yleft + (yright - yleft) * (10. * std::pow(normalized_time, 3.) - 15. * std::pow(normalized_time, 4) +
                                       6. * std::pow(normalized_time, 5.));
}

ScalarType smoothstep_c3(ScalarType normalized_time, ScalarType yleft, ScalarType yright)
{
    return yleft + (yright - yleft) * (35. * std::pow(normalized_time, 4.) - 84. * std::pow(normalized_time, 5.) +
                                       70. * std::pow(normalized_time, 6.) - 20. * std::pow(normalized_time, 7.));
}

ScalarType smoothstep_c4(ScalarType normalized_time, ScalarType yleft, ScalarType yright)
{
    return yleft + (yright - yleft) * (126. * std::pow(normalized_time, 5.) - 420. * std::pow(normalized_time, 6.) +
                                       540. * std::pow(normalized_time, 7.) - 315. * std::pow(normalized_time, 8.) +
                                       70. * std::pow(normalized_time, 9.));
}

std::vector<ScalarType> linspace(ScalarType xleft, ScalarType xright, size_t num_steps)
{
    std::vector<ScalarType> result(num_steps);
    if (num_steps == 1) {
        result[0] = xleft;
        return result;
    }
    const ScalarType step = (xright - xleft) / static_cast<ScalarType>(num_steps - 1);
    for (size_t i = 0; i < num_steps; ++i) {
        result[i] = xleft + static_cast<ScalarType>(i) * step;
    }
    return result;
}

std::vector<ScalarType>
evaluate_function(const std::vector<ScalarType>& evaluation_points,
                  const std::function<ScalarType(ScalarType, ScalarType, ScalarType)>& smoother_function)
{
    using namespace params;

    std::vector<ScalarType> results(evaluation_points.size());

    for (size_t i = 0; i < evaluation_points.size(); ++i) {
        ScalarType normalized_time = (evaluation_points[i] - xleft) / (xright - xleft);
        results[i]                 = smoother_function(normalized_time, yleft, yright);
    }

    return results;
}

mio::IOResult<void> run_benchmarks(std::string result_dir, size_t num_runs, size_t num_warm_up_runs)
{
    using namespace params;

    std::vector<ScalarType> smoothercos_times(num_runs);
    std::vector<ScalarType> smoothstep_c1_times(num_runs);
    std::vector<ScalarType> smoothstep_c2_times(num_runs);
    std::vector<ScalarType> smoothstep_c3_times(num_runs);
    std::vector<ScalarType> smoothstep_c4_times(num_runs);
    std::vector<ScalarType> evaluation_times = linspace(xleft, xright, num_steps);

    const auto benchmark = [&](const std::function<ScalarType(ScalarType, ScalarType, ScalarType)>& smoother_function,
                               std::vector<ScalarType>& times) {
        for (size_t run = 0; run < num_warm_up_runs; ++run) {
            const auto results = evaluate_function(evaluation_times, smoother_function);
            (void)results;
        }

        for (size_t run = 0; run < num_runs; ++run) {
            mio::timing::BasicTimer timer;
            timer.start();
            const auto results = evaluate_function(evaluation_times, smoother_function);
            timer.stop();
            (void)results;
            times[run] = mio::timing::time_in_seconds(timer.get_elapsed_time());
        }
    };

    benchmark(smoothercos, smoothercos_times);
    benchmark(smoothstep_c1, smoothstep_c1_times);
    benchmark(smoothstep_c2, smoothstep_c2_times);
    benchmark(smoothstep_c3, smoothstep_c3_times);
    benchmark(smoothstep_c4, smoothstep_c4_times);

    mio::TimeSeries<ScalarType> timeseries(5);
    for (size_t i = 0; i < num_runs; ++i) {
        Eigen::VectorXd values(5);
        values << smoothercos_times[i], smoothstep_c1_times[i], smoothstep_c2_times[i], smoothstep_c3_times[i],
            smoothstep_c4_times[i];
        timeseries.add_time_point(i, values);
    }

    std::string filename = fmt::format("smoother_times_{}warmupruns_{}runs_t0={}_tmax={}_numsteps={}.csv",
                                       num_warm_up_runs, num_runs, xleft, xright, num_steps);
    BOOST_OUTCOME_TRY(
        timeseries.export_csv(mio::path_join(result_dir, filename),
                              {"smoothercos", "smoothstep_c1", "smoothstep_c2", "smoothstep_c3", "smoothstep_c4"}));

    return mio::success();
}

int main(int argc, char** argv)
{

    auto cli_parameters =
        mio::cli::ParameterSetBuilder()
            .add<"ResultDirectory">(mio::path_join(std::filesystem::current_path().string(), "results_runtime"))
            .add<"NumberRuns">(100, {.alias = "nRun"})
            .add<"NumberWarmupRuns">(20, {.alias = "nWURun"})
            .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, cli_parameters, {"ResultDirectory"});
    if (!cli_result) {
        std::cout << cli_result.error().message();
        return cli_result.error().code().value();
    }

    std::filesystem::path res_dir(cli_parameters.get<"ResultDirectory">());
    std::filesystem::create_directories(res_dir);

    auto result = run_benchmarks(cli_parameters.get<"ResultDirectory">(), cli_parameters.get<"NumberRuns">(),
                                 cli_parameters.get<"NumberWarmupRuns">());
    return 0;
}