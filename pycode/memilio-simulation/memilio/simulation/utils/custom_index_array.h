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

#include "pybind11/detail/common.h"
#include "pybind_util.h"
#include "utils/index.h"
#include "memilio/utils/index.h"
#include "memilio/utils/custom_index_array.h"

#include "pybind11/pybind11.h"
#include "pybind11/pytypes.h"
#include "pybind11/cast.h"

#include "boost/optional.hpp"
#include <type_traits>
#include <stdexcept>

namespace pymio
{

// Recursively bind the members of custom index array
// that require a single Tag as a template argument.
template <class C>
void bind_single_tag_template_members(pybind11::module_& m, pybind11::class_<C>&)
{
}
template <class C, class T, class... Ts>
void bind_single_tag_template_members(pybind11::module_& m, pybind11::class_<C>& c)
{
    std::string tname = pretty_name<T>();
    c.def(("size_" + tname).c_str(), &C::template size<T>);

    // Catch warning ImportError: generic_type: type "" is already registered!
    try {
        bind_Index<T>(m, ("Index_" + tname).c_str());
    }
    catch (std::runtime_error& e) {
    }
    bind_single_tag_template_members<C, Ts...>(m, c); //next Tag
}

// Represents a python index slice, i.e. start:stop:step
// where each argument can be None.
template <class Index>
struct SliceExpression {
    boost::optional<Index> start;
    boost::optional<Index> stop;
    boost::optional<Index> step;
};

// Convert any index expression into a slice expression.
// The expression is either a python slice expression (i.e. `start:stop:step`)
// or a single index.
template <class Tag>
SliceExpression<mio::Index<Tag>> convert_to_slice_expr(const pybind11::object& obj)
{
    if (pybind11::isinstance<pybind11::slice>(obj)) {
        //convert a python slice expression
        auto slice = obj.cast<pybind11::slice>();
        return {cast_or_none<mio::Index<Tag>>(slice.attr("start")), cast_or_none<mio::Index<Tag>>(slice.attr("stop")),
                cast_or_none<mio::Index<Tag>>(slice.attr("step"))};
    }
    else {
        //make a slice from a single index, i.e. `i:i+1:1`
        auto idx = obj.cast<mio::Index<Tag>>();
        return {idx, idx + mio::Index<Tag>{1}, mio::Index<Tag>{1}};
    }
}

// Recursively slices an array along each dimension.
// C: original array
// S: slice of previous dimensions
// T / Ts: category tags of the remaining dimensions of the array
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
    const auto I          = C::Index::size - sizeof...(Ts) - 1;
    const auto slice_expr = convert_to_slice_expr<T>(tup[I]);
    const auto size       = self.template size<T>();

    const auto start = get_optional_value_or(slice_expr.start, mio::Index<T>(0));
    if (start >= size) {
        //invalid start is the only part that may cause an exception
        throw pybind11::index_error("Out of range.");
    }
    //python slices just ignore everything beyond the end of the array, but C++ crashes, so constrain stop in [start, size]
    const auto stop = std::min(std::max(get_optional_value_or(slice_expr.stop, size), start), size);
    const auto step = get_optional_value_or(slice_expr.step, mio::Index<T>(1));
    const auto n    = (size_t(stop) - size_t(start)) / size_t(step); //C++ uses number of elements instead of end index

    auto array_slice = slice.template slice<T>({size_t(start), n, size_t(step)});
    return get_slice_of_array<C, decltype(array_slice), Ts...>(self, array_slice, tup); //recurse to next dimension
}

// Assign scalar to array.
// self: custom index array
// indices: python slice expression or single index
// value: scalar value to assign
template <class C, class T, class... Tags>
void assign_scalar(C& self, const pybind11::object& indices, const T& value)
{
    const auto index_tuple = self.size().size == 1 ? pybind11::make_tuple(indices) : indices.cast<pybind11::tuple>();
    if (index_tuple.size() != self.size().size) {
        throw pybind11::index_error("Invalid number of dimensions.");
    }
    get_slice_of_array<C, C, Tags...>(self, self, index_tuple) = value;
}

// Bind members that are different for arrays with single or multi index.
template <class C, class Tag>
void bind_single_or_multi_index_members_CustomIndexArray(pybind11::module_& m, pybind11::class_<C>& c,
                                                         std::string const& name)
{
    c.def("size", [](const C& self) {
        return self.size(); //just a single index, no tuple
    });
}
template <class C, class... Tags>
std::enable_if_t<(sizeof...(Tags) > 1)> bind_single_or_multi_index_members_CustomIndexArray(pybind11::module_& m,
                                                                                            pybind11::class_<C>& c,
                                                                                            std::string const& name)
{
    c.def("size", [](const C& self) {
        return self.size().indices; //tuple of single indices
    });

    // Catch warning ImportError: generic_type: type "" is already registered!
    try {
        bind_MultiIndex<Tags...>(m, ("MultiIndex_" + name).c_str());
    }
    catch (std::runtime_error& e) {
    }
}

template <class Type, class... Tags>
void bind_CustomIndexArray(pybind11::module_& m, std::string const& name)
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
        .def("get_flat_index", &C::get_flat_index)
        //scalar assignment
        .def("__setitem__", &assign_scalar<C, typename C::value_type, Tags...>)
        //scalar assignment with conversion from double
        //TODO: may need SFINAE in the future, only compiles if value type is convertible from double, e.g. UncertainValue
        .def("__setitem__", &assign_scalar<C, double, Tags...>);
    //TODO: __setitem__ with list or numpy array, e.g. array[AgeGroup(0):AgeGroup(3)] = [1, 2, 3]
    //TODO: __getitem__. Is it ever necessary to store a reference to a slice?

    bind_single_or_multi_index_members_CustomIndexArray<C, Tags...>(m, c, name);
    bind_single_tag_template_members<C, Tags...>(m, c);
}

} // namespace pymio

#endif //PYMIO_CUSTOM_INDEX_ARRAY_H
