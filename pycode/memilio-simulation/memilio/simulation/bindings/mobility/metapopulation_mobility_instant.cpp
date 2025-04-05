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
#include "mobility/metapopulation_mobility_instant.h"
#include "pybind_util.h"

#include "pybind11/eigen.h"

namespace py = pybind11;

namespace pymio
{

void bind_mobility_parameters(py::module_& m, std::string const& name)
{
    bind_class<mio::MobilityParameters<double>, EnablePickling::IfAvailable>(m, name.c_str())
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def(py::init<const mio::MobilityCoefficientGroup&>(), py::arg("coeffs"))
        .def_property(
            "coefficients", py::overload_cast<>(&mio::MobilityParameters<double>::get_coefficients),
            [](mio::MobilityParameters<double>& self, const mio::MobilityCoefficientGroup& v) {
                self.get_coefficients() = v;
            },
            py::return_value_policy::reference_internal);
}

void bind_mobility_parameter_edge(py::module_& m, std::string const& name)
{
    bind_class<mio::Edge<mio::MobilityParameters<double>>, EnablePickling::IfAvailable>(m, name.c_str())
        .def_property_readonly("start_node_idx",
                               [](const mio::Edge<mio::MobilityParameters<double>>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const mio::Edge<mio::MobilityParameters<double>>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property",
            [](const mio::Edge<mio::MobilityEdge<double>>& self) -> auto& {
                return self.property;
            },
            py::return_value_policy::reference_internal);
}

void bind_mobility(py::module_& m, std::string const& name)
{
    bind_class<mio::MobilityEdge<double>, EnablePickling::IfAvailable>(m, name.c_str())
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def(py::init<const mio::MobilityParameters<double>&>(), py::arg("params"))
        .def_property_readonly(
            "parameters",
            [](const mio::MobilityEdge<double>& self) -> auto& {
                return self.get_parameters();
            },
            py::return_value_policy::reference_internal);
}

void bind_mobility_edge(py::module_& m, std::string const& name)
{
    bind_class<mio::Edge<mio::MobilityEdge<double>>, EnablePickling::IfAvailable>(m, name.c_str())
        .def_property_readonly("start_node_idx",
                               [](const mio::Edge<mio::MobilityEdge<double>>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const mio::Edge<mio::MobilityEdge<double>>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property",
            [](const mio::Edge<mio::MobilityEdge<double>>& self) -> auto& {
                return self.property;
            },
            py::return_value_policy::reference_internal);
}

} // namespace pymio
