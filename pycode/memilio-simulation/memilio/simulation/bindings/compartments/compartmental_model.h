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
#ifndef PYMIO_COMPARTMENTALMODEL_H
#define PYMIO_COMPARTMENTALMODEL_H

#include "memilio/compartments/compartmental_model.h"
#include "pybind_util.h"

#include "pybind11/pybind11.h"

namespace pymio
{

/*
 * @brief bind a compartmental model for any Populations and Parameters class
 */
template <class InfectionState, class Populations, class Parameters, EnablePickling F>
void bind_CompartmentalModel(pybind11::module_& m, std::string const& name)
{
    using Model = mio::CompartmentalModel<double, InfectionState, Populations, Parameters>;
    bind_class<Model, F>(m, name.c_str())
        .def(pybind11::init<Populations const&, Parameters const&>())
        .def("apply_constraints", &Model::apply_constraints)
        .def("check_constraints", &Model::check_constraints)
        .def("get_initial_values", &Model::get_initial_values)
        .def_property(
            "populations",
            [](const Model& self) -> auto& {
                return self.populations;
            },
            [](Model& self, Populations& p) {
                self.populations = p;
            })
        .def_property(
            "parameters",
            [](const Model& self) -> auto& {
                return self.parameters;
            },
            [](Model& self, Parameters& p) {
                self.parameters = p;
            });
}

} // namespace pymio

#endif //PYMIO_COMPARTMENTALMODEL_H
