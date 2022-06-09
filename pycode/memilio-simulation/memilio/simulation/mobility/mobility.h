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
#ifndef PYMIO_MOBILITY_H
#define PYMIO_MOBILITY_H

#include "memilio/mobility/mobility.h"
#include "memilio/mobility/graph.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{

template <class Simulation>
void bind_MigrationGraph(py::module& m, std::string const& name)
{
    using G = mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>;
    py::class_<G>(m, name.c_str())
        .def(py::init<>())
        .def(
            "add_node", [](G & self, int id, const typename Simulation::Model& p, double t0, double dt) -> auto& {
                return self.add_node(id, p, t0, dt);
            },
            py::arg("id"), py::arg("model"), py::arg("t0") = 0.0, py::arg("dt") = 0.1,
            py::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const mio::MigrationParameters&>, py::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const Eigen::VectorXd&>, py::return_value_policy::reference_internal)
        .def_property_readonly("num_nodes",
                               [](const G& self) {
                                   return self.nodes().size();
                               })
        .def(
            "get_node", [](const G& self, size_t node_idx) -> auto& { return self.nodes()[node_idx]; },
            py::return_value_policy::reference_internal)
        .def_property_readonly("num_edges",
                               [](const G& self) {
                                   return self.edges().size();
                               })
        .def(
            "get_edge", [](const G& self, size_t edge_idx) -> auto& { return self.edges()[edge_idx]; },
            py::return_value_policy::reference_internal)
        .def("get_num_out_edges",
             [](const G& self, size_t node_idx) {
                 return self.out_edges(node_idx).size();
             })
        .def(
            "get_out_edge",
            [](const G& self, size_t node_idx, size_t edge_idx) -> auto& { return self.out_edges(node_idx)[edge_idx]; },
            py::return_value_policy::reference_internal);
}

void bind_migration_parameters(py::module& m, std::string const& name);

void bind_migration_parameter_edge(py::module& m, std::string const& name);

void bind_migration(py::module& m, std::string const& name);

void bind_migration_edge(py::module& m, std::string const& name);

template <typename Model>
void bind_ModelNode(py::module& m, std::string const& name)
{
    py::class_<mio::Node<Model>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const mio::Node<Model>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property", [](const mio::Node<Model>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);
}

template <typename Simulation>
void bind_SimulationNode(py::module& m, std::string const& name)
{
    py::class_<mio::Node<mio::SimulationNode<Simulation>>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const mio::Node<Simulation>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property",
            [](const mio::Node<mio::SimulationNode<Simulation>>& self) -> auto& { return self.property.get_simulation(); },
            py::return_value_policy::reference_internal);
}

/*
 * @brief bind Graph for any node and edge type
 */
template <class Model>
void bind_ModelGraph(py::module& m, std::string const& name)
{
    using G = mio::Graph<Model, mio::MigrationParameters>;
    py::class_<G>(m, name.c_str())
        .def(py::init<>())
        .def("add_node", &G::template add_node<const Model&>, py::arg("id"), py::arg("model"), py::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const mio::MigrationParameters&>, py::arg("start_node_idx"), py::arg("end_node_idx"), py::arg("migration_parameters"),
             py::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const Eigen::VectorXd&>, py::return_value_policy::reference_internal)
        .def_property_readonly("num_nodes",
                               [](const G& self) {
                                   return self.nodes().size();
                               })
        .def(
            "get_node", [](const G& self, size_t node_idx) -> auto& { return self.nodes()[node_idx]; },
            py::return_value_policy::reference_internal)
        .def_property_readonly("num_edges",
                               [](const G& self) {
                                   return self.edges().size();
                               })
        .def(
            "get_edge", [](const G& self, size_t edge_idx) -> auto& { return self.edges()[edge_idx]; },
            py::return_value_policy::reference_internal)
        .def("get_num_out_edges",
             [](const G& self, size_t node_idx) {
                 return self.out_edges(node_idx).size();
             })
        .def(
            "get_out_edge",
            [](const G& self, size_t node_idx, size_t edge_idx) -> auto& { return self.out_edges(node_idx)[edge_idx]; },
            py::return_value_policy::reference_internal);
}

} // namespace pymio

#endif //PYMIO_MOBILITY_H