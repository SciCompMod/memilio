/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Maximilian Betz
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
#include "mobility/mobility.h"
#include "memilio/mobility/graph.h"

#include "pybind11/eigen.h"

namespace pymio
{

void bind_migration_parameters(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::MigrationParameters>(m, name.c_str())
        .def(pybind11::init<const Eigen::VectorXd&>(), pybind11::arg("coeffs"))
        .def(pybind11::init<const mio::MigrationCoefficientGroup&>(), pybind11::arg("coeffs"))
        .def_property(
            "coefficients", pybind11::overload_cast<>(&mio::MigrationParameters::get_coefficients),
            [](mio::MigrationParameters& self, const mio::MigrationCoefficientGroup& v) {
                self.get_coefficients() = v;
            },
            pybind11::return_value_policy::reference_internal);
}

void bind_migration_parameter_edge(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::Edge<mio::MigrationParameters>>(m, name.c_str())
        .def_property_readonly("start_node_idx",
                               [](const mio::Edge<mio::MigrationParameters>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const mio::Edge<mio::MigrationParameters>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property", [](const mio::Edge<mio::MigrationEdge>& self) -> auto& { return self.property; },
            pybind11::return_value_policy::reference_internal);
}

void bind_migration(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::MigrationEdge>(m, name.c_str())
        .def(pybind11::init<const Eigen::VectorXd&>(), pybind11::arg("coeffs"))
        .def(pybind11::init<const mio::MigrationParameters&>(), pybind11::arg("params"))
        .def_property_readonly(
            "parameters", [](const mio::MigrationEdge& self) -> auto& { return self.get_parameters(); },
            pybind11::return_value_policy::reference_internal);
}

void bind_migration_edge(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::Edge<mio::MigrationEdge>>(m, name.c_str())
        .def_property_readonly("start_node_idx",
                               [](const mio::Edge<mio::MigrationEdge>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const mio::Edge<mio::MigrationEdge>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property", [](const mio::Edge<mio::MigrationEdge>& self) -> auto& { return self.property; },
            pybind11::return_value_policy::reference_internal);
}

} // namespace pymio