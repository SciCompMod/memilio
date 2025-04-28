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
#include "utils/uncertain_value.h"
#include "pybind_util.h"
#include "memilio/utils/uncertain_value.h"

namespace py = pybind11;

namespace pymio
{

void bind_uncertain_value(py::module_& m, std::string const& name)
{
    bind_class<mio::UncertainValue<double>, EnablePickling::Required>(m, name.c_str())
        .def(py::init<double>(), py::arg("value") = 0.0)
        .def_property(
            "value",
            [](mio::UncertainValue<double>& self) {
                return ScalarType(self);
            },
            [](mio::UncertainValue<double>& self, double v) {
                self = v;
            })
        .def("set_distribution", //a property would be nicer but getter and setter use a different type
             &mio::UncertainValue<double>::set_distribution)
        .def(
            "get_distribution",
            [](const mio::UncertainValue<double>& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def(
            "get_distribution",
            [](mio::UncertainValue<double>& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def("draw_sample", &mio::UncertainValue<double>::draw_sample);
}

} // namespace pymio
