#ifndef PYMIO_PARAMETER_SET_H
#define PYMIO_PARAMETER_SET_H

#include "memilio/utils/parameter_set.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{

template <class ParameterSet>
auto bind_ParameterSet(py::module& m, std::string const& name)
{
    py::class_<ParameterSet> c(m, name.c_str());
    mio::foreach_tag<ParameterSet>([&c](auto t) {
        using Tag = decltype(t);

        //CAUTION: This requires ParameterTag::name() to be unique within the ParameterSet
        c.def_property(
            Tag::name().c_str(), [](const ParameterSet& self) -> auto& { return self.template get<Tag>(); },
            [](ParameterSet& self, typename Tag::Type const& v) {
                self.template get<Tag>() = v;
            },
            py::return_value_policy::reference_internal);
    });
    return c;
}

} // namespace pymio

#endif //PYMIO_PARAMETER_SET_H