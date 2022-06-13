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
#ifndef PYMIO_CUSTOM_INDEX_ARRAY_H
#define PYMIO_CUSTOM_INDEX_ARRAY_H

#include "pybind_util.h"
#include "memilio/utils/index.h"
#include "memilio/utils/custom_index_array.h"

#include "pybind11/pybind11.h"

namespace pymio
{

template <class C>
void bind_templated_members_CustomIndexArray(pybind11::class_<C>&)
{
}

template <class C, class T, class... Ts>
void bind_templated_members_CustomIndexArray(pybind11::class_<C>& c)
{
    std::string tname = pretty_name<T>();
    c.def(("size_" + tname).c_str(), &C::template size<T>);

    // recursively bind the member for each type
    bind_templated_members_CustomIndexArray<C, Ts...>(c);
}

template <class Type, class... Tags>
void bind_CustomIndexArray(pybind11::module& m, std::string const& name)
{
    using C     = typename mio::CustomIndexArray<Type, Tags...>;
    using Index = typename mio::CustomIndexArray<Type, Tags...>::Index;
    pybind11::class_<C> c(m, name.c_str());
    c.def(pybind11::init([](Index const& sizes, Type const& val) {
         return C(sizes, val);
     }))
        .def(pybind11::init([](Index const& sizes) {
            return C(sizes);
        }))
        .def("numel", &C::numel)
        .def(
            "__getitem__", [](const C& self, Index const& idx) -> auto& { return self[idx]; },
            pybind11::return_value_policy::reference_internal)
        .def(
            "__getitem__",
            [](const C& self,
               std::tuple<mio::Index<Tags>...> idx) -> auto& { //python natively handles multi-indices as tuples
                return self[{std::get<mio::Index<Tags>>(idx)...}];
            },
            pybind11::return_value_policy::reference_internal)
        .def("__setitem__",
             [](C& self, Index const& idx, double value) {
                 self[idx] = value;
             })
        .def("__setitem__",
             [](C& self, std::tuple<mio::Index<Tags>...> idx, double value) {
                 self[{std::get<mio::Index<Tags>>(idx)...}] = value;
             })
        .def(
            "__iter__",
            [](const C& s) {
                return pybind11::make_iterator(s.begin(), s.end());
            },
            pybind11::keep_alive<0, 1>())
        .def("get_flat_index", &C::get_flat_index);

    // Not supported in Python yet: Slicing

    // bind all templated members for types in Tags...
    bind_templated_members_CustomIndexArray<C, Tags...>(c);
}

} // namespace pymio

#endif //PYMIO_CUSTOM_INDEX_ARRAY_H