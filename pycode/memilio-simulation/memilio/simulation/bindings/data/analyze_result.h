/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Martin Siggel, Daniel Abele, Martin J. Kuehn, Jan Kleinert, Maximilian Betz
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
#ifndef PYMIO_DATA_ANALYZE_RESULT_H
#define PYMIO_DATA_ANALYZE_RESULT_H

#include "memilio/data/analyze_result.h"
#include "pybind_util.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{

/**
 * @brief Bind functions interpolate_simulation_result and interpolate_ensemble_results.
 */
void bind_interpolate_result_methods(py::module_& m)
{

    m.def(
        "interpolate_simulation_result",
        [](const mio::TimeSeries<double>& ts) {
            return mio::interpolate_simulation_result(ts);
        },
        py::arg("ts"));
    m.def("interpolate_simulation_result",
          static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&, const double)>(
              &mio::interpolate_simulation_result),
          py::arg("ts"), py::arg("abs_tol"));

    m.def("interpolate_simulation_result",
          static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&, const std::vector<double>&)>(
              &mio::interpolate_simulation_result),
          py::arg("ts"), py::arg("interpolation_times"));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<mio::TimeSeries<double>>);
}

} // namespace pymio

#endif // PYMIO_DATA_ANALYZE_RESULT_H
