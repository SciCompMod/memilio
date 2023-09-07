/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Henrik Zunker
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
#ifndef PYMIO_EPI_DYNAMIC_LOCKDOWN_H
#define PYMIO_EPI_DYNAMIC_LOCKDOWN_H

#include "memilio/epidemiology/dynamic_npis.h"
#include "pybind_util.h"

#include "pybind11/pybind11.h"

namespace pymio
{

void bind_dynamicNPI_members(pybind11::module_& m, std::string const& name)
{
    using C    = mio::DynamicNPIs;
    py::class_<C>(m, name.c_str())
        .def(py::init<>())
        .def("set_interval", &C::set_interval)
        .def("set_duration", &C::set_duration)
        .def("set_base_value", &C::set_base_value)
        .def("set_threshold", &C::set_threshold)
}

} // namespace pymio

#endif //PYMIO_EPI_DYNAMIC_LOCKDOWN_H
