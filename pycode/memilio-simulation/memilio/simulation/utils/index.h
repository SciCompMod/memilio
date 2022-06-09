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