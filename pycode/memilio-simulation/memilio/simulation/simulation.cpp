/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin Siggel, Daniel Abele, Martin J. Kuehn, Jan Kleinert
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

#include "pybind_util.h"
#include "epidemiology/damping.h"
#include "epidemiology/contact_matrix.h"
#include "epidemiology/damping_sampling.h"
#include "utils/date.h"
#include "utils/time_series.h"
#include "utils/parameter_distributions.h"
#include "utils/uncertain_value.h"
#include "memilio/mobility/mobility.h"
#include "memilio/utils/date.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/regions.h"
#include "memilio/epidemiology/uncertain_matrix.h"

namespace py = pybind11;

PYBIND11_MODULE(_simulation, m)
{
    pymio::bind_date(m, "Date");

    auto damping_class = py::class_<mio::SquareDamping>(m, "Damping");
    pymio::bind_damping_members(damping_class);

    auto dampings_class = py::class_<mio::SquareDampings>(m, "Dampings");
    pymio::bind_dampings_members(dampings_class);
    
    pymio::bind_time_series(m, "TimeSeries");

    pymio::bind_parameter_distribution(m, "ParameterDistribution");
    pymio::bind_parameter_distribution_normal(m, "ParameterDistributionNormal");
    pymio::bind_parameter_distribution_uniform(m, "ParameterDistributionUniform");

    pymio::bind_uncertain_value(m, "UncertainValue");

    auto contact_matrix_class = py::class_<mio::ContactMatrix>(m, "ContactMatrix");
    pymio::bind_damping_expression_members(contact_matrix_class);
    contact_matrix_class.def_property_readonly("num_groups", &mio::ContactMatrix::get_num_groups);

    auto contact_matrix_group_class = py::class_<mio::ContactMatrixGroup>(m, "ContactMatrixGroup");
    pymio::bind_damping_expression_group_members(contact_matrix_group_class);
    contact_matrix_group_class.def_property_readonly("num_groups", &mio::ContactMatrixGroup::get_num_groups);

    pymio::bind_damping_sampling(m, "DampingSampling");

    py::class_<mio::UncertainContactMatrix>(m, "UncertainContactMatrix")
        .def(py::init<>())
        .def(py::init<const mio::ContactMatrixGroup&>())
        .def_property(
            "cont_freq_mat", py::overload_cast<>(&mio::UncertainContactMatrix::get_cont_freq_mat),
            [](mio::UncertainContactMatrix& self, const mio::ContactMatrixGroup& c) {
                self.get_cont_freq_mat() = c;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "dampings", py::overload_cast<>(&mio::UncertainContactMatrix::get_dampings),
            [](mio::UncertainContactMatrix& self, const std::vector<mio::DampingSampling>& v) {
                self.get_dampings() = v;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "school_holidays",
            [](const mio::UncertainContactMatrix& self) {
                std::vector<std::pair<double, double>> v(self.get_school_holidays().size());
                std::transform(self.get_school_holidays().begin(), self.get_school_holidays().end(), v.begin(),
                               [](auto& p) {
                                   return std::make_pair(double(p.first), double(p.second));
                               });
                return v;
            },
            [](mio::UncertainContactMatrix& self, const std::vector<std::pair<double, double>>& v) {
                self.get_school_holidays().resize(v.size());
                std::transform(v.begin(), v.end(), self.get_school_holidays().begin(), [](auto& p) {
                    return std::make_pair(mio::SimulationTime(p.first), mio::SimulationTime(p.second));
                });
            })
        .def_property("school_holiday_damping",
                      py::overload_cast<>(&mio::UncertainContactMatrix::get_school_holiday_damping),
                      [](mio::UncertainContactMatrix& self, const mio::DampingSampling& v) {
                          self.get_school_holiday_damping() = v;
                      });

    auto migration_damping_class = py::class_<mio::VectorDamping>(m, "MigrationDamping");
    pymio::bind_damping_members(migration_damping_class);

    auto migration_dampings_class = py::class_<mio::VectorDampings>(m, "MigrationDampings");
    pymio::bind_dampings_members(migration_dampings_class);

    auto migration_coeffs_class = py::class_<mio::MigrationCoefficients>(m, "MigrationCoefficients");
    pymio::bind_damping_expression_members(migration_coeffs_class);

    auto migration_coeff_group_class = py::class_<mio::MigrationCoefficientGroup>(m, "MigrationCoefficientGroup");
    pymio::bind_damping_expression_group_members(migration_coeff_group_class);

    py::class_<mio::MigrationParameters>(m, "MigrationParameters")
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def(py::init<const mio::MigrationCoefficientGroup&>(), py::arg("coeffs"))
        .def_property(
            "coefficients", py::overload_cast<>(&mio::MigrationParameters::get_coefficients),
            [](mio::MigrationParameters& self, const mio::MigrationCoefficientGroup& v) {
                self.get_coefficients() = v;
            },
            py::return_value_policy::reference_internal);

    py::class_<mio::Edge<mio::MigrationParameters>>(m, "MigrationParameterEdge")
        .def_property_readonly("start_node_idx",
                               [](const mio::Edge<mio::MigrationParameters>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const mio::Edge<mio::MigrationParameters>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property", [](const mio::Edge<mio::MigrationEdge>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);

    py::class_<mio::MigrationEdge>(m, "Migration")
        .def(py::init<const Eigen::VectorXd&>(), py::arg("coeffs"))
        .def(py::init<const mio::MigrationParameters&>(), py::arg("params"))
        .def_property_readonly(
            "parameters", [](const mio::MigrationEdge& self) -> auto& { return self.get_parameters(); },
            py::return_value_policy::reference_internal);

    py::class_<mio::Edge<mio::MigrationEdge>>(m, "MigrationEdge")
        .def_property_readonly("start_node_idx",
                               [](const mio::Edge<mio::MigrationEdge>& self) {
                                   return self.start_node_idx;
                               })
        .def_property_readonly("end_node_idx",
                               [](const mio::Edge<mio::MigrationEdge>& self) {
                                   return self.end_node_idx;
                               })
        .def_property_readonly(
            "property", [](const mio::Edge<mio::MigrationEdge>& self) -> auto& { return self.property; },
            py::return_value_policy::reference_internal);

    m.def(
        "get_state_id_de",
        [](int county) {
            return int(mio::regions::de::get_state_id(mio::regions::de::CountyId(county)));
        },
        py::arg("county_id"));
    m.def(
        "get_holidays_de",
        [](int state, mio::Date start_date, mio::Date end_date) {
            auto h = mio::regions::de::get_holidays(mio::regions::de::StateId(state), start_date, end_date);
            return std::vector<std::pair<mio::Date, mio::Date>>(h.begin(), h.end());
        },
        py::arg("state_id"), py::arg("start_date") = mio::Date(std::numeric_limits<int>::min(), 1, 1),
        py::arg("end_date") = mio::Date(std::numeric_limits<int>::max(), 1, 1));

    py::enum_<mio::LogLevel>(m, "LogLevel")
        .value("Off", mio::LogLevel::off)
        .value("Critical", mio::LogLevel::critical)
        .value("Error", mio::LogLevel::err)
        .value("Warning", mio::LogLevel::warn)
        .value("Info", mio::LogLevel::info)
        .value("Debug", mio::LogLevel::debug)
        .value("Trace", mio::LogLevel::trace);
    m.def("set_log_level", &mio::set_log_level);

    m.attr("__version__") = "dev";
}
