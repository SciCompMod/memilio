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
#include "epidemiology/damping_sampling.h"
#include "pybind_util.h"
#include "memilio/epidemiology/damping_sampling.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/utils/uncertain_value.h"

namespace pymio
{

void bind_damping_sampling(pybind11::module& m, std::string const& name)
{
    pymio::pybind_pickle_class<mio::DampingSampling>(m, name.c_str())
        .def(pybind11::init([](const mio::UncertainValue& value, int level, int type, double time,
                         const std::vector<size_t>& matrices, const Eigen::Ref<const Eigen::VectorXd>& groups) {
                 return mio::DampingSampling(value, mio::DampingLevel(level), mio::DampingType(type),
                                             mio::SimulationTime(time), matrices, groups);
             }),
             pybind11::arg("value"), pybind11::arg("level"), pybind11::arg("type"), pybind11::arg("time"), pybind11::arg("matrix_indices"),
             pybind11::arg("group_weights"))
        .def_property("value", pybind11::overload_cast<>(&mio::DampingSampling::get_value), &mio::DampingSampling::set_value,
                      pybind11::return_value_policy::reference_internal)
        .def_property(
            "level",
            [](const mio::DampingSampling& self) {
                return int(self.get_level());
            },
            [](mio::DampingSampling& self, int lvl) {
                self.set_level(mio::DampingLevel(lvl));
            })
        .def_property(
            "type",
            [](const mio::DampingSampling& self) {
                return int(self.get_type());
            },
            [](mio::DampingSampling& self, int typ) {
                self.set_type(mio::DampingType(typ));
            })
        .def_property(
            "time",
            [](const mio::DampingSampling& self) {
                return double(self.get_time());
            },
            [](mio::DampingSampling& self, double t) {
                self.set_time(mio::SimulationTime(t));
            })
        .def_property("matrix_indices", &mio::DampingSampling::get_matrix_indices,
                      &mio::DampingSampling::set_matrix_indices)
        .def_property(
            "group_weights", &mio::DampingSampling::get_group_weights,
            [](mio::DampingSampling& self, const Eigen::Ref<const Eigen::VectorXd>& v) {
                self.set_group_weights(v);
            },
            pybind11::return_value_policy::reference_internal);
}

} // namespace pymio