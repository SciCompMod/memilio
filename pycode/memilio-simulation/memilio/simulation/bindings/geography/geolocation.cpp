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
#ifndef PYMIO_GEOLOCATION_H
#define PYMIO_GEOLOCATION_H

#include "memilio/math/integrator.h"
#include "memilio/geography/geolocation.h"
#include "pybind_util.h"

#include "pybind11/pybind11.h"
#include <pybind11/eigen.h>

namespace py = pybind11;

namespace pymio
{

void bind_geolocation(pybind11::module_& m, std::string const& name)
{
    bind_class<mio::geo::GeographicalLocation, EnablePickling::Never>(m, name.c_str()).def(py::init<double, double>());
}
} // namespace pymio

#endif //PYMIO_GEOLOCATION_H