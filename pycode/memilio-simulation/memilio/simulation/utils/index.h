#ifndef PYMIO_INDEX_H
#define PYMIO_INDEX_H

#include "memilio/utils/index.h"

#include "pybind11/pybind11.h"
#include "pybind11/operators.h"

namespace pymio
{

// bind an index for a single tag
template <class Tag> 
void bind_Index(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::Index<Tag>> c(m, name.c_str());
    c.def(pybind11::init<size_t>(), pybind11::arg("value"));
    c.def(pybind11::self == pybind11::self);
    c.def(pybind11::self != pybind11:: self);
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
void bind_MultiIndex(pybind11::module& m, std::string const& name)
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