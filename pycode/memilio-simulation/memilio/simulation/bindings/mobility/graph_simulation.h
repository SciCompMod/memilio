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
#ifndef PYMIO_GRAPH_SIMULATION_H
#define PYMIO_GRAPH_SIMULATION_H

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/graph_simulation.h"
#include "pybind_util.h"

#include "pybind11/pybind11.h"

namespace pymio
{

/*
 * @brief bind GraphSimulation for any node and edge type
 */
template <class Graph>
void bind_GraphSimulation(pybind11::module_& m, std::string const& name)
{
    using GS = mio::GraphSimulation<Graph>;
    bind_class<GS, EnablePickling::Never>(m, name.c_str())
        .def(pybind11::init([](Graph& graph, double t0, double dt) {
                 return std::make_unique<GS>(mio::make_mobility_sim(t0, dt, std::move(graph)));
             }),
             pybind11::arg("graph"), pybind11::arg("t0") = 0.0, pybind11::arg("dt") = 1.0)
        .def_property_readonly(
            "graph",
            [](GS& self) -> Graph& {
                return self.get_graph();
            },
            pybind11::return_value_policy::reference_internal)
        .def_property_readonly("t", &GS::get_t)
        .def("advance", &GS::advance, pybind11::arg("tmax"));
}

} // namespace pymio

#endif //PYMIO_GRAPH_SIMULATION_H
