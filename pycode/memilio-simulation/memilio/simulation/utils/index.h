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
#ifndef PYMIO_INDEX_H
#define PYMIO_INDEX_H

#include "memilio/utils/index.h"

#include "pybind11/pybind11.h"
#include "pybind11/operators.h"

namespace py = pybind11;

namespace pymio
{

// bind an index for a single tag
template <class Tag> 
void bind_Index(py::module& m, std::string const& name)
{
    py::class_<mio::Index<Tag>> c(m, name.c_str());
    c.def(py::init<size_t>(), py::arg("value"));
    c.def(py::self == py::self);
    c.def(py::self != py:: self);
}

// helper function for implicitly casting from py::tuple to Index in Python.
// This extracts an Index from a py::tuple of Indices from the correct position,
// given the corresponding Index type
template <typename Tag, class Tuple>
mio::Index<Tag> extract_index(py::tuple& t)
{
    return t[mio::details::IndexPosition<Tag, Tuple>::value].template cast<mio::Index<Tag>>();
}

// bind an index for more than one tag
template <class... Tags>
void bind_MultiIndex(py::module& m, std::string const& name)
{
    using C = mio::Index<Tags...>;
    py::class_<C> c(m, name.c_str());
    c.def(py::init<mio::Index<Tags> const&...>()).def(py::init([](py::tuple t) {
        return C(extract_index<Tags, C>(t)...);
    }));

    py::implicitly_convertible<py::tuple, C>();
}


} // namespace pymio

#endif //PYMIO_INDEX_H