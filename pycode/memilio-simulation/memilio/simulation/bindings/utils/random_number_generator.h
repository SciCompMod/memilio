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

#include "memilio/utils/random_number_generator.h"
#include "pybind_util.h"

namespace py = pybind11;

namespace pymio
{
void bind_random_number_generator(py::module_& m, std::string const& name)
{
    bind_class<mio::RandomNumberGenerator, EnablePickling::Never>(m, name.c_str())
        .def(py::init<>())
        .def_property_readonly("key", &mio::RandomNumberGenerator::get_key)
        .def_property("counter", &mio::RandomNumberGenerator::get_counter, &mio::RandomNumberGenerator::set_counter)
        .def_property_readonly("seeds", &mio::RandomNumberGenerator::get_seeds)
        .def("increment_counter", &mio::RandomNumberGenerator::increment_counter)
        .def("seed", &mio::RandomNumberGenerator::seed, py::arg("seeds"));
}

void bind_discrete_distribution(py::module_& m, std::string const& name)
{
    bind_class<mio::DiscreteDistribution<int>, EnablePickling::Never>(m, name.c_str())
        .def_static("get_instance", &mio::DiscreteDistribution<int>::get_instance,
                    py::return_value_policy::
                        reference) // Note: reference_internal cannot be used here because this is a static function
        .def("__call__", [](mio::DiscreteDistribution<int>& self, mio::RandomNumberGenerator& rng,
                            std::vector<double> distribution) {
            return self(rng, distribution);
        });
}
} // namespace pymio