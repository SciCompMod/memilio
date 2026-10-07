/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Kilian Volmer, Henrik Zunker
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
#include "pybind_util.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "abm_halle_pais.cpp"

namespace py = pybind11;


PYBIND11_MODULE(_simulation_halle_abm, m)
{
    m.def(
    "run_once",
    [](std::string person_file, std::string contact_dir, std::string history_file,
        int history_lookback_days, bool allow_missing_history, const std::vector<double>& theta,
        mio::Date start_date, int num_days) {
        auto result = run_once(person_file, contact_dir, history_file, history_lookback_days,
                                allow_missing_history, theta, start_date, num_days);
        if (!result) {
            throw std::runtime_error(result.error().formatted_message());
        }
        const auto& observable = result.value();
        return py::make_tuple(observable.days, observable.values);
    },
    "Run the ABM once. Returns (days, values) with one [deaths, pais_medium, pais_severe] entry per day.",
    py::arg("person_file"), py::arg("contact_dir"), py::arg("history_file"),
    py::arg("history_lookback_days"), py::arg("allow_missing_history"), py::arg("theta"),
    py::arg("start_date"), py::arg("num_days"));

    m.attr("__version__") = "dev";
}
