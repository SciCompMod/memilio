/* 
* Copyright (C) 2020-2026 MEmilio
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
#include "utils/date.h"
#include "pybind_util.h"
#include "memilio/utils/date.h"

namespace py = pybind11;

namespace pymio
{

void bind_date(py::module_& m, std::string const& name)
{
    bind_class<mio::Date, EnablePickling::Required>(m, name.c_str())
        .def(py::init<int, int, int>(), py::arg("year"), py::arg("month"), py::arg("day"))
        .def_readwrite("year", &mio::Date::year)
        .def_readwrite("month", &mio::Date::month)
        .def_readwrite("day", &mio::Date::day)
        .def_property_readonly("day_in_year",
                               [](const mio::Date& self) {
                                   return mio::get_day_in_year(self);
                               })
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def(py::self <= py::self)
        .def(py::self > py::self)
        .def(py::self >= py::self)
        // Add offset to a date (new_date = date + offset)
        .def("__add__",
             [](const mio::Date& self, int offset_days) {
                 return mio::offset_date_by_days(self, offset_days);
             })
        // In-place add offset to a date (date += offset)
        .def("__iadd__",
             [](mio::Date& self, int offset_days) {
                 self = mio::offset_date_by_days(self, offset_days);
                 return self;
             })
        // Substract offset from a date (new_date = date - offset)
        .def("__sub__",
             [](const mio::Date& self, int offset_days) {
                 return mio::offset_date_by_days(self, -offset_days);
             })
        // In-place substract offset from a date (date -= offset)
        .def("__isub__",
             [](mio::Date& self, int offset_days) {
                 self = mio::offset_date_by_days(self, -offset_days);
                 return self;
             })
        // Substract dare from date to get their difference (difference = date1 - date2)
        .def("__sub__",
             [](const mio::Date& self, const mio::Date& other) {
                 return mio::get_offset_in_days(self, other);
             })
        .def("__str__", [](const mio::Date& self) {
            return std::to_string(self.year) + "." + std::to_string(self.month) + "." + std::to_string(self.day);
        });
}

} // namespace pymio
