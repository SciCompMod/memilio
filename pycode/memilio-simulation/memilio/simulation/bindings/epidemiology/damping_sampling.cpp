/* 
* Copyright (C) 2020-2026 MEmilio
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
#include "epidemiology/damping_sampling.h"
#include "pybind_util.h"
#include "memilio/epidemiology/damping_sampling.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/utils/uncertain_value.h"

#include "pybind11/eigen.h"

namespace py = pybind11;

namespace pymio
{

void bind_damping_sampling(py::module_& m, std::string const& name)
{
    bind_class<mio::DampingSampling<double>, EnablePickling::Required>(m, name.c_str())
        .def(py::init([](const mio::UncertainValue<double>& value, int level, int type, double time,
                         const std::vector<size_t>& matrices, const Eigen::Ref<const Eigen::VectorXd>& groups) {
                 return mio::DampingSampling(value, mio::DampingLevel(level), mio::DampingType(type),
                                             mio::SimulationTime<double>(time), matrices, groups);
             }),
             py::arg("value"), py::arg("level"), py::arg("type"), py::arg("time"), py::arg("matrix_indices"),
             py::arg("group_weights"))
        .def_property("value", py::overload_cast<>(&mio::DampingSampling<double>::get_value),
                      &mio::DampingSampling<double>::set_value, py::return_value_policy::reference_internal)
        .def_property(
            "level",
            [](const mio::DampingSampling<double>& self) {
                return int(self.get_level());
            },
            [](mio::DampingSampling<double>& self, int lvl) {
                self.set_level(mio::DampingLevel(lvl));
            })
        .def_property(
            "type",
            [](const mio::DampingSampling<double>& self) {
                return int(self.get_type());
            },
            [](mio::DampingSampling<double>& self, int typ) {
                self.set_type(mio::DampingType(typ));
            })
        .def_property(
            "time",
            [](const mio::DampingSampling<double>& self) {
                return double(self.get_time());
            },
            [](mio::DampingSampling<double>& self, double t) {
                self.set_time(mio::SimulationTime<double>(t));
            })
        .def_property("matrix_indices", &mio::DampingSampling<double>::get_matrix_indices,
                      &mio::DampingSampling<double>::set_matrix_indices)
        .def_property(
            "group_weights", &mio::DampingSampling<double>::get_group_weights,
            [](mio::DampingSampling<double>& self, const Eigen::Ref<const Eigen::VectorXd>& v) {
                self.set_group_weights(v);
            },
            py::return_value_policy::reference_internal);
}

} // namespace pymio
