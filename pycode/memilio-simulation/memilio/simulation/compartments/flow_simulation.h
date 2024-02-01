/* 
* Copyright (C) 2020-2024 MEmilio
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
#ifndef PYMIO_FLOW_SIMULATION_H
#define PYMIO_FLOW_SIMULATION_H

#include "memilio/compartments/flow_simulation.h"

#include "pybind11/pybind11.h"

namespace pymio
{

template <class Model>
void bind_Flow_Simulation(pybind11::module_& m)
{
    pybind11::class_<mio::FlowSimulation<Model>>(m, "FlowSimulation")
        .def(pybind11::init<const Model&, double, double>(), pybind11::arg("model"), pybind11::arg("t0") = 0,
             pybind11::arg("dt") = 0.1)
        .def_property_readonly("result",
                               pybind11::overload_cast<>(&mio::FlowSimulation<Model>::get_result, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def_property_readonly("flows", &mio::FlowSimulation<Model>::get_flows,
                               pybind11::return_value_policy::reference_internal)
        .def("advance", &mio::FlowSimulation<Model>::advance, pybind11::arg("tmax"));
}

} // namespace pymio

#endif //PYMIO_FLOW_SIMULATION_H
