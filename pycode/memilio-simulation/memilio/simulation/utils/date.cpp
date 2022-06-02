#include "utils/date.h"

namespace pymio
{

void bind_date(pybind11::module& m, std::string const& name)
{
    pymio::pybind_pickle_class<mio::Date>(m, name.c_str())
        .def(pybind11::init<int, int, int>(), pybind11::arg("year"), pybind11::arg("month"), pybind11::arg("day"))
        .def_readwrite("year", &mio::Date::year)
        .def_readwrite("month", &mio::Date::month)
        .def_readwrite("day", &mio::Date::day)
        .def(pybind11::self == pybind11::self)
        .def(pybind11::self != pybind11::self);
}

} // namespace pymio