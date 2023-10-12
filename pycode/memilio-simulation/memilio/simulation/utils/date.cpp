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
#include "utils/date.h"
#include "pybind_util.h"
#include "memilio/utils/date.h"

namespace py = pybind11;

namespace pymio
{

void bind_date(py::module_& m, std::string const& name)
{
    pymio::pybind_pickle_class<mio::Date>(m, name.c_str())
        .def(py::init<int, int, int>(), py::arg("year"), py::arg("month"), py::arg("day"))
        .def_readwrite("year", &mio::Date::year)
        .def_readwrite("month", &mio::Date::month)
        .def_readwrite("day", &mio::Date::day)
        .def(py::self == py::self)
        .def(py::self != py::self);
}

} // namespace pymio
