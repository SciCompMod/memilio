/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Carlotta Gerstein
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

#include "memilio/math/time_series_functor.h"
#include "memilio/math/interpolation.h"

#include "pybind_util.h"
#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{

void bind_time_series_functor(py::module_& m, std::string const& name)
{
    bind_class<mio::TimeSeriesFunctor<double>, EnablePickling::Never>(m, name.c_str())
        .def(py::init())
        .def(py::init<mio::TimeSeriesFunctorType, mio::TimeSeries<double>>())
        .def(py::init([](const mio::TimeSeries<double>& data) {
            return mio::TimeSeriesFunctor(mio::TimeSeriesFunctorType::LinearInterpolation, data);
        }))
        .def(py::init([](std::vector<std::vector<double>>&& table) {
            return mio::TimeSeriesFunctor<double>(mio::TimeSeriesFunctorType::LinearInterpolation, table);
        }))
        .def("__call__", [](mio::TimeSeriesFunctor<double>& self, double time) {
            return self(time);
        });
}

} // namespace pymio
