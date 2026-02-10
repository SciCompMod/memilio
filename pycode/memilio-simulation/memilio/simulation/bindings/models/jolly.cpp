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

// Macro to skip main() in jolly_bindings.cpp
#define JOLLY_BINDINGS_SKIP_MAIN
#include "jolly_bindings.cpp"

namespace py = pybind11;

PYBIND11_MODULE(_simulation_jolly, m)
{
    m.def("simulate", &simulate, "Simulates the jolly model", py::arg("farm_file"), py::arg("tmax") = 40,
          py::arg("dt") = 1.0, py::arg("suspicion_threshold") = 0.05, py::arg("sensitivity") = 0.95,
          py::arg("h0") = 0.001, py::arg("r0") = 1000, py::arg("alpha") = 2, py::arg("A0_SEI") = 0.2,
          py::arg("A0_EI") = 0.2, py::arg("A0_ID") = 0.2, py::arg("A0_DeathRate") = 0.000128, py::arg("A1_SEI") = 0.2,
          py::arg("A1_EI") = 0.2, py::arg("A1_ID") = 0.2, py::arg("A1_DeathRate") = 0.000128, py::arg("A2_SEI") = 0.2,
          py::arg("A2_EI") = 0.2, py::arg("A2_ID") = 0.2, py::arg("A2_DeathRate") = 0.000588, py::arg("A3_SEI") = 0.2,
          py::arg("A3_EI") = 0.2, py::arg("A3_ID") = 0.2, py::arg("A3_DeathRate") = 0.000588, py::arg("A4_SEI") = 0.2,
          py::arg("A4_EI") = 0.2, py::arg("A4_ID") = 0.2, py::arg("A4_DeathRate") = 0.000588,
          py::arg("foi_inner_factor0") = 1.0, py::arg("foi_outer_factor0") = 1.0, py::arg("foi_inner_factor1") = 1.0,
          py::arg("foi_outer_factor1") = 1.0, py::arg("foi_inner_factor2") = 1.0, py::arg("foi_outer_factor2") = 1.0,
          py::arg("foi_inner_factor3") = 1.0, py::arg("foi_outer_factor3") = 1.0, py::arg("foi_inner_factor4") = 1.0,
          py::arg("foi_outer_factor4") = 1.0, py::arg("first_infection_day") = 0, py::arg("second_infection_day") = 2,
          py::arg("third_infection_day") = 2, py::arg("seed") = 42);

    m.attr("__version__") = "dev";
}
