#ifndef PYMIO_INDEX_H
#define PYMIO_INDEX_H

#include "memilio/utils/index.h"

#include "pybind11/pybind11.h"

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

} // namespace pymio

#endif //PYMIO_INDEX_H