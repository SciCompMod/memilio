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
#include "memilio/utils/custom_index_array.h"

#include "pybind11/pybind11.h"
#include "pybind11/operators.h"

namespace pymio
{

//bind members of Index class that only exist if the Tag is not an enum type
template <class Tag>
std::enable_if_t<!std::is_enum<Tag>::value> bind_Index_members_if_enum(pybind11::class_<mio::Index<Tag>>& c)
{
}

//bind members of Index class that only exist if the Tag is an enum type
template <class Tag>
std::enable_if_t<std::is_enum<Tag>::value> bind_Index_members_if_enum(pybind11::class_<mio::Index<Tag>>& c)
{
    pybind11::implicitly_convertible<Tag, mio::Index<Tag>>();
}

// bind an index for a single tag
template <class Tag>
void bind_Index(pybind11::module_& m, std::string const& name)
{
    pybind11::class_<mio::Index<Tag>> c(m, name.c_str());
    c.def(pybind11::init<size_t>(), pybind11::arg("value"));
    c.def(pybind11::self == pybind11::self);
    c.def(pybind11::self != pybind11::self);

    bind_Index_members_if_enum(c);
}

// helper function for implicitly casting from pybind11::tuple to Index in Python.
// This extracts an Index from a pybind11::tuple of Indices from the correct position,
// given the corresponding Index type
template <typename Tag, class Tuple>
mio::Index<Tag> extract_index(pybind11::tuple& t)
{
    return t[mio::details::IndexPosition<Tag, Tuple>::value].template cast<mio::Index<Tag>>();
}

// bind an index for more than one tag
template <class... Tags>
void bind_MultiIndex(pybind11::module_& m, std::string const& name)
{
    using C = mio::Index<Tags...>;
    pybind11::class_<C> c(m, name.c_str());
    c.def(pybind11::init<mio::Index<Tags> const&...>()).def(pybind11::init([](pybind11::tuple t) {
        return C(extract_index<Tags, C>(t)...);
    }));

    pybind11::implicitly_convertible<pybind11::tuple, C>();
}

} // namespace pymio

#endif //PYMIO_INDEX_H
