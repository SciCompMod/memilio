/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef PYMIO_PARAMETER_DISTRIBUTIONS_H
#define PYMIO_PARAMETER_DISTRIBUTIONS_H

#include "pybind11/pybind11.h"

namespace pymio
{

void bind_parameter_distribution(pybind11::module_& m, std::string const& name);

void bind_parameter_distribution_normal(pybind11::module_& m, std::string const& name);

void bind_parameter_distribution_uniform(pybind11::module_& m, std::string const& name);

} // namespace pymio

#endif //PYMIO_PARAMETER_DISTRIBUTIONS_H
