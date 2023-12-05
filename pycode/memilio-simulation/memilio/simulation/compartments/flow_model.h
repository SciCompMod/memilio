/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: 
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
#ifndef PYMIO_FLOW_MODEL_H
#define PYMIO_FLOW_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/compartments/flow_model.h"

#include "pybind11/pybind11.h"
#include <string>

namespace pymio
{

template <class Flows, class Model, class BindType>
void bind_FlowIndex(BindType& c)
{
    mio::foreach_tag<Flows>([&c](auto t) {
        using Tag = decltype(t);

        c.def_property_readonly(std::to_string(Tag::source) + "_" + std::to_string(Tag::target),
                                [](const Model& self, const typename Model::FlowIndex& index) -> size_t {
                                    return self.get_flat_flow_index<Tag::source, Tag::target>(index);
                                });
    });
}

/*
 * @brief bind a Flow model for any Populations and Parameters class
 */
template <class InfectionState, class Populations, class Parameters, class Flows>
void bind_FlowModel(pybind11::module_& m, std::string const& name)
{
    using Model = mio::FlowModel<InfectionState, Populations, Parameters, Flows>;
    pybind11::class_<Model, mio::CompartmentalModel<InfectionState, Populations, Parameters>> c(m, name.c_str());
    c.def(pybind11::init<Populations const&, Parameters const&>());

    bind_FlowIndex<Flows, Model>(c);
}

} // namespace pymio

#endif //PYMIO_FLOW_MODEL_H
