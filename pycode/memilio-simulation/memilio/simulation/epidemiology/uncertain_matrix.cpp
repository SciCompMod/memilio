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
#include "epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/damping_sampling.h"

#include <pybind11/stl.h>

namespace py = pybind11;

namespace pymio
{

void bind_uncertain_contact_matrix(py::module_& m, std::string const& name)
{
    py::class_<mio::UncertainContactMatrix>(m, name.c_str())
        .def(py::init<>())
        .def(py::init<const mio::ContactMatrixGroup&>())
        .def_property(
            "cont_freq_mat", py::overload_cast<>(&mio::UncertainContactMatrix::get_cont_freq_mat),
            [](mio::UncertainContactMatrix& self, const mio::ContactMatrixGroup& c) {
                self.get_cont_freq_mat() = c;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "dampings", py::overload_cast<>(&mio::UncertainContactMatrix::get_dampings),
            [](mio::UncertainContactMatrix& self, const std::vector<mio::DampingSampling>& v) {
                self.get_dampings() = v;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "school_holidays",
            [](const mio::UncertainContactMatrix& self) {
                std::vector<std::pair<double, double>> v(self.get_school_holidays().size());
                std::transform(self.get_school_holidays().begin(), self.get_school_holidays().end(), v.begin(),
                               [](auto& p) {
                                   return std::make_pair(double(p.first), double(p.second));
                               });
                return v;
            },
            [](mio::UncertainContactMatrix& self, const std::vector<std::pair<double, double>>& v) {
                self.get_school_holidays().resize(v.size());
                std::transform(v.begin(), v.end(), self.get_school_holidays().begin(), [](auto& p) {
                    return std::make_pair(mio::SimulationTime(p.first), mio::SimulationTime(p.second));
                });
            })
        .def_property("school_holiday_damping",
                      py::overload_cast<>(&mio::UncertainContactMatrix::get_school_holiday_damping),
                      [](mio::UncertainContactMatrix& self, const mio::DampingSampling& v) {
                          self.get_school_holiday_damping() = v;
                      });
}

} // namespace pymio
