/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "pybind_util.h"

#include "pybind11/pybind11.h"

namespace pymio
{

/*
 * @brief bind Simulation for any number model
 */
template <class Simulation>
void bind_Simulation(pybind11::module_& m, std::string const& name)
{

    bind_class<Simulation, EnablePickling::IfAvailable>(m, name.c_str())
        .def(pybind11::init<const typename Simulation::Model&, double, double>(), pybind11::arg("model"),
             pybind11::arg("t0") = 0, pybind11::arg("dt") = 0.1)
        .def_property_readonly("result", pybind11::overload_cast<>(&Simulation::get_result, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def_property_readonly("model", pybind11::overload_cast<>(&Simulation::get_model, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def_property_readonly("dt", pybind11::overload_cast<>(&Simulation::get_dt, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def_property("integrator", pybind11::overload_cast<>(&Simulation::get_integrator, pybind11::const_),
                               &Simulation::set_integrator, pybind11::return_value_policy::reference_internal)
        .def("advance", &Simulation::advance, pybind11::arg("tmax"))
        .doc() = "A class for the simulation of a compartment model.";
}

} // namespace pymio

#endif //PYMIO_SIMULATION_H
