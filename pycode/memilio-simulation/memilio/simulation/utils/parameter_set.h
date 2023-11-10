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
#ifndef PYMIO_PARAMETER_SET_H
#define PYMIO_PARAMETER_SET_H

#include "memilio/utils/parameter_set.h"

#include "pybind11/pybind11.h"

namespace pymio
{

template <class ParameterSet>
auto bind_ParameterSet(pybind11::module_& m, std::string const& name)
{
    pybind11::class_<ParameterSet> c(m, name.c_str());
    mio::foreach_tag<ParameterSet>([&c](auto t) {
        using Tag = decltype(t);

        //CAUTION: This requires ParameterTag::name() to be unique within the ParameterSet
        c.def_property(
            Tag::name().c_str(), [](const ParameterSet& self) -> auto& { return self.template get<Tag>(); },
            [](ParameterSet& self, typename Tag::Type const& v) {
                self.template get<Tag>() = v;
            },
            pybind11::return_value_policy::reference_internal);
    });
    return c;
}

} // namespace pymio

#endif //PYMIO_PARAMETER_SET_H
