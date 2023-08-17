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

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/graph.h"

#include "pybind11/pybind11.h"

namespace pymio
{

template <class Simulation>
void bind_MigrationGraph(pybind11::module_& m, std::string const& name)
{
    using G = mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>;
    pybind11::class_<G>(m, name.c_str())
        .def(pybind11::init<>())
        .def(
            "add_node",
            [](G& self, int id, const typename Simulation::Model& p, double t0, double dt) -> auto& {
                return self.add_node(id, p, t0, dt);
            },
            pybind11::arg("id"), pybind11::arg("model"), pybind11::arg("t0") = 0.0, pybind11::arg("dt") = 0.1,
            pybind11::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const mio::MigrationParameters&>,
             pybind11::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const Eigen::VectorXd&>,
             pybind11::return_value_policy::reference_internal)
        .def_property_readonly("num_nodes",
                               [](const G& self) {
                                   return self.nodes().size();
                               })
        .def(
            "get_node",
            [](const G& self, size_t node_idx) -> auto& {
                return self.nodes()[node_idx];
            },
            pybind11::return_value_policy::reference_internal)
        .def_property_readonly("num_edges",
                               [](const G& self) {
                                   return self.edges().size();
                               })
        .def(
            "get_edge",
            [](const G& self, size_t edge_idx) -> auto& {
                return self.edges()[edge_idx];
            },
            pybind11::return_value_policy::reference_internal)
        .def("get_num_out_edges",
             [](const G& self, size_t node_idx) {
                 return self.out_edges(node_idx).size();
             })
        .def(
            "get_out_edge",
            [](const G& self, size_t node_idx, size_t edge_idx) -> auto& {
                return self.out_edges(node_idx)[edge_idx];
            },
            pybind11::return_value_policy::reference_internal);
}

void bind_migration_parameters(pybind11::module_& m, std::string const& name);

void bind_migration_parameter_edge(pybind11::module_& m, std::string const& name);

void bind_migration(pybind11::module_& m, std::string const& name);

void bind_migration_edge(pybind11::module_& m, std::string const& name);

template <typename Model>
void bind_ModelNode(pybind11::module_& m, std::string const& name)
{
    pybind11::class_<mio::Node<Model>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const mio::Node<Model>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property",
            [](const mio::Node<Model>& self) -> auto& {
                return self.property;
            },
            pybind11::return_value_policy::reference_internal);
}

template <typename Simulation>
void bind_SimulationNode(pybind11::module_& m, std::string const& name)
{
    pybind11::class_<mio::Node<mio::SimulationNode<Simulation>>>(m, name.c_str())
        .def_property_readonly("id",
                               [](const mio::Node<Simulation>& self) {
                                   return self.id;
                               })
        .def_property_readonly(
            "property",
            [](const mio::Node<mio::SimulationNode<Simulation>>& self) -> auto& {
                return self.property.get_simulation();
            },
            pybind11::return_value_policy::reference_internal);
}

/*
 * @brief bind Graph for any node and edge type
 */
template <class Model>
void bind_ModelGraph(pybind11::module_& m, std::string const& name)
{
    using G = mio::Graph<Model, mio::MigrationParameters>;
    pybind11::class_<G>(m, name.c_str())
        .def(pybind11::init<>())
        .def("add_node", &G::template add_node<const Model&>, pybind11::arg("id"), pybind11::arg("model"),
             pybind11::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const mio::MigrationParameters&>, pybind11::arg("start_node_idx"),
             pybind11::arg("end_node_idx"), pybind11::arg("migration_parameters"),
             pybind11::return_value_policy::reference_internal)
        .def("add_edge", &G::template add_edge<const Eigen::VectorXd&>,
             pybind11::return_value_policy::reference_internal)
        .def_property_readonly("num_nodes",
                               [](const G& self) {
                                   return self.nodes().size();
                               })
        .def(
            "get_node",
            [](const G& self, size_t node_idx) -> auto& {
                return self.nodes()[node_idx];
            },
            pybind11::return_value_policy::reference_internal)
        .def_property_readonly("num_edges",
                               [](const G& self) {
                                   return self.edges().size();
                               })
        .def(
            "get_edge",
            [](const G& self, size_t edge_idx) -> auto& {
                return self.edges()[edge_idx];
            },
            pybind11::return_value_policy::reference_internal)
        .def("get_num_out_edges",
             [](const G& self, size_t node_idx) {
                 return self.out_edges(node_idx).size();
             })
        .def(
            "get_out_edge",
            [](const G& self, size_t node_idx, size_t edge_idx) -> auto& {
                return self.out_edges(node_idx)[edge_idx];
            },
            pybind11::return_value_policy::reference_internal);
}

} // namespace pymio

#endif //PYMIO_MOBILITY_H
