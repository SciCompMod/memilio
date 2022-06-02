/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Maximilian Betz
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
#include "utils/parameter_distributions.h"
#include "pybind_util.h"
#include "memilio/utils/parameter_distributions.h"

namespace pymio
{

void bind_parameter_distribution(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::ParameterDistribution>(m, name.c_str())
        .def_property("lower_bound", &mio::ParameterDistribution::get_lower_bound,
                      &mio::ParameterDistribution::set_lower_bound)
        .def_property("upper_bound", &mio::ParameterDistribution::get_upper_bound,
                      &mio::ParameterDistribution::set_upper_bound)
        .def("add_predefined_sample", &mio::ParameterDistribution::add_predefined_sample)
        .def("remove_predefined_samples", &mio::ParameterDistribution::remove_predefined_samples)
        .def("get_sample", &mio::ParameterDistribution::get_sample);
}

void bind_parameter_distribution_normal(pybind11::module& m, std::string const& name)
{
    pymio::pybind_pickle_class<mio::ParameterDistributionNormal, mio::ParameterDistribution>(m, name.c_str())
        .def(pybind11::init<double, double, double, double>(), pybind11::arg("lb"), pybind11::arg("ub"), pybind11::arg("mean"),
             pybind11::arg("std_dev"))
        .def(pybind11::init<double, double, double>(), pybind11::arg("lb"), pybind11::arg("ub"), pybind11::arg("mean"))
        .def_property("mean", &mio::ParameterDistributionNormal::get_mean, &mio::ParameterDistributionNormal::set_mean)
        .def_property("standard_dev", &mio::ParameterDistributionNormal::get_standard_dev,
                      &mio::ParameterDistributionNormal::set_standard_dev);
}

void bind_parameter_distribution_uniform(pybind11::module& m, std::string const& name)
{
    pymio::pybind_pickle_class<mio::ParameterDistributionUniform, mio::ParameterDistribution>(m, name.c_str())
        .def(pybind11::init<>())
        .def(pybind11::init<double, double>(), pybind11::arg("lb"), pybind11::arg("ub"));
}

} // namespace pymio