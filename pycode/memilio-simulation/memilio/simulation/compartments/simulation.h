/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#ifndef PYMIO_SIMULATION_H
#define PYMIO_SIMULATION_H

#include "memilio/compartments/simulation.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{

/*
 * @brief bind Simulation for any number model
 */
template <class Simulation>
void bind_Simulation(py::module& m, std::string const& name)
{
    py::class_<Simulation>(m, name.c_str())
        .def(py::init<const typename Simulation::Model&, double, double>(), py::arg("model"), py::arg("t0") = 0,
             py::arg("dt") = 0.1)
        .def_property_readonly("result", py::overload_cast<>(&Simulation::get_result, py::const_),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("model", py::overload_cast<>(&Simulation::get_model, py::const_),
                               py::return_value_policy::reference_internal)
        .def("advance", &Simulation::advance, py::arg("tmax"));
}

} // namespace pymio

#endif //PYMIO_SIMULATION_H