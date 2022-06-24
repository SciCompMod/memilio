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
#include "pybind11/pytypes.h"
#include "pybind11/cast.h"

#include "boost/optional.hpp"

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

template<class T, class Obj>
boost::optional<T> cast_if_not_none(Obj&& obj)
{
    if (obj.is_none()) {
        return {};
    } else {
        return obj.template cast<T>();
    }
}

template<class Index>
struct Slice
{
    boost::optional<Index> start;
    boost::optional<Index> stop;
    boost::optional<Index> step;
};

//convert an index expression into a slice.
//the expression is either a python slice expression (i.e. `start:stop:step`)
//or a single index.
template <class Tag>
Slice<mio::Index<Tag>> make_slice(const pybind11::object& obj)
{
    if (pybind11::isinstance<pybind11::slice>(obj)) {
        //convert a python slice expression
        auto slice = obj.cast<pybind11::slice>();
        return {cast_if_not_none<mio::Index<Tag>>(slice.attr("start")),
                cast_if_not_none<mio::Index<Tag>>(slice.attr("stop")),
                cast_if_not_none<mio::Index<Tag>>(slice.attr("step"))};
    }
    else {
        //make a slice from a single index, i.e. `i:i+1:1`
        auto idx = obj.cast<mio::Index<Tag>>();
        return {idx, idx + mio::Index<Tag>{1}, mio::Index<Tag>{1}};
    }
}

//recursively slices an array along each dimension
//I: index of current dimension
//C: original array
//S: slice of previous dimensions
//T / Ts: category tags of the remaining dimensions of the array
template <class C, class S>
auto get_slice_of_array(C& /*self*/, S& slice, const pybind11::tuple& /*tup*/)
{
    //base case of recursion, no more dimensions, slice finished
    return slice;
}
template <class C, class S, class T, class... Ts>
auto get_slice_of_array(C& self, S& slice, const pybind11::tuple& tup)
{
    //normal case of recursion
    //slice along dimension I, identified by tag T
    const auto I = C::Index::size - sizeof...(Ts) - 1;
    auto c_slice = make_slice<T>(tup[I]);
    auto size    = mio::get<I>(self.size());
    auto start   = get_optional_value_or(c_slice.start, mio::Index<T>(0));
    auto stop    = get_optional_value_or(c_slice.stop, size);
    auto step    = get_optional_value_or(c_slice.step, mio::Index<T>(1));
    if (start >= size) {
        throw pybind11::index_error("Out of range.");
    }
    stop   = std::min(std::max(stop, start), size); //python slicing just ignores everything beyond the end of the array
    auto n = (size_t(stop) - size_t(start)) / size_t(step);
    auto array_slice = slice.template slice<T>({size_t(start), n, size_t(step)});
    return get_slice_of_array<C, decltype(array_slice), Ts...>(self, array_slice, tup);
}

//bind slicing for one dimensional arrays
//accepts single slice as indices
template<class C, class... Ts>
void bind_slicing_operations_CustomIndexArray(pybind11::class_<C>& c) 
{
    static_assert(sizeof...(Ts) == C::Index::size, "");
    c.def("__setitem__", [](C& self, pybind11::object indices, const typename C::value_type& value) {
        auto index_tuple = self.size().size == 1 ? pybind11::make_tuple(indices) : indices.cast<pybind11::tuple>();
        if (index_tuple.size() != self.size().size) {
            throw pybind11::index_error("Invalid number of dimensions.");
        }
        get_slice_of_array<C, C, Ts...>(self, self, index_tuple) = value;
    });
    //TODO: set slice from (numpy) array
    //TODO: __getitem__
}

// //bind slicing for one dimensional arrays
// //accepts tuple of slices as indices
// template<class C, class T, class... Ts>
// std::enable_if_t<(sizeof...(Ts) > 0)> bind_slicing_operations_CustomIndexArray(pybind11::class_<C>& c) 
// {    
//     c.def("__setitem__", [](C& self, pybind11::tuple tup, const typename C::value_type& value) {
//         if (tup.size() != self.size().size) {
//             throw pybind11::index_error("Invalid number of dimensions.");
//         }
//         get_slice_of_array<0, C, C, T, Ts...>(self, self, tup) = value;
//     });
// }

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

    // slicing of arrays
    bind_slicing_operations_CustomIndexArray<C, Tags...>(c);

    // bind all members of CustomIndexArray that work on a single parameter
    bind_templated_members_CustomIndexArray<C, Tags...>(c);
}

} // namespace pymio

#endif //PYMIO_CUSTOM_INDEX_ARRAY_H