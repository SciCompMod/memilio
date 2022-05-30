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

#include "templates.h"
#include "pybind_util.h"
#include "epidemiology/damping.h"
#include "epidemiology/contact_matrix.h"
#include "memilio/utils/date.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/regions.h"
#include "memilio/epidemiology/uncertain_matrix.h"

namespace py = pybind11;

PYBIND11_MODULE(_simulation, m)
{
    pymio::pybind_pickle_class<mio::Date>(m, "Date")
        .def(py::init<int, int, int>(), py::arg("year"), py::arg("month"), py::arg("day"))
        .def_readwrite("year", &mio::Date::year)
        .def_readwrite("month", &mio::Date::month)
        .def_readwrite("day", &mio::Date::day)
        .def(py::self == py::self)
        .def(py::self != py::self);

    auto damping_class = py::class_<mio::SquareDamping>(m, "Damping");
    pymio::bind_damping_members(damping_class);

    auto dampings_class = py::class_<mio::SquareDampings>(m, "Dampings");
    pymio::bind_dampings_members(dampings_class);

    py::class_<mio::TimeSeries<double>>(m, "TimeSeries")
        .def(py::init<Eigen::Index>(), py::arg("num_elements"))
        .def("get_num_time_points", &mio::TimeSeries<double>::get_num_time_points)
        .def("get_num_elements", &mio::TimeSeries<double>::get_num_elements)
        .def("get_time", py::overload_cast<Eigen::Index>(&mio::TimeSeries<double>::get_time), py::arg("index"))
        .def("get_last_time", py::overload_cast<>(&mio::TimeSeries<double>::get_last_time))
        .def("get_value", py::overload_cast<Eigen::Index>(&mio::TimeSeries<double>::get_value), py::arg("index"))
        .def("get_last_value", py::overload_cast<>(&mio::TimeSeries<double>::get_last_value))
        .def("__len__", &mio::TimeSeries<double>::get_num_time_points)
        .def(
            "__getitem__",
            [](mio::TimeSeries<double>& self, Eigen::Index i) {
                if (i >= 0 && i < self.get_num_time_points()) {
                    return self[i];
                }
                else {
                    throw pybind11::index_error("Index out of range."); //needs to throw exception for iterable
                }
            },
            py::is_operator(), py::arg("index"))
        .def(
            "__setitem__",
            [](mio::TimeSeries<double>& self, Eigen::Index i, Eigen::Ref<const mio::TimeSeries<double>::Vector> expr) {
                if (i >= 0 && i < self.get_num_time_points()) {
                    self[i] = expr;
                }
                else {
                    throw pybind11::index_error("Index out of range."); //needs to throw exception for iterable
                }
            },
            py::is_operator(), py::arg("index"), py::arg("v"))
        .def("add_time_point",
             [](mio::TimeSeries<double>& self) {
                 return self.add_time_point();
             })
        .def("add_time_point",
             [](mio::TimeSeries<double>& self, double t) {
                 return self.add_time_point(t);
             })
        .def("add_time_point",
             [](mio::TimeSeries<double>& self, double t, Eigen::Ref<const mio::TimeSeries<double>::Vector> expr) {
                 return self.add_time_point(t, expr);
             })
        .def("as_ndarray", [](mio::TimeSeries<double>& self) {
            auto m = Eigen::Map<mio::TimeSeries<double>::Matrix>(self.data(), self.get_num_rows(),
                                                                 self.get_num_time_points());
            return Eigen::Ref<mio::TimeSeries<double>::Matrix>(m);
        });

    py::class_<mio::ParameterDistribution>(m, "ParameterDistribution")
        .def_property("lower_bound", &mio::ParameterDistribution::get_lower_bound,
                      &mio::ParameterDistribution::set_lower_bound)
        .def_property("upper_bound", &mio::ParameterDistribution::get_upper_bound,
                      &mio::ParameterDistribution::set_upper_bound)
        .def("add_predefined_sample", &mio::ParameterDistribution::add_predefined_sample)
        .def("remove_predefined_samples", &mio::ParameterDistribution::remove_predefined_samples)
        .def("get_sample", &mio::ParameterDistribution::get_sample);

    pymio::pybind_pickle_class<mio::ParameterDistributionNormal, mio::ParameterDistribution>(m, "ParameterDistributionNormal")
        .def(py::init<double, double, double, double>(), py::arg("lb"), py::arg("ub"), py::arg("mean"),
             py::arg("std_dev"))
        .def(py::init<double, double, double>(), py::arg("lb"), py::arg("ub"), py::arg("mean"))
        .def_property("mean", &mio::ParameterDistributionNormal::get_mean, &mio::ParameterDistributionNormal::set_mean)
        .def_property("standard_dev", &mio::ParameterDistributionNormal::get_standard_dev,
                      &mio::ParameterDistributionNormal::set_standard_dev);

    pymio::pybind_pickle_class<mio::ParameterDistributionUniform, mio::ParameterDistribution>(m, "ParameterDistributionUniform")
        .def(py::init<>())
        .def(py::init<double, double>(), py::arg("lb"), py::arg("ub"));

    pymio::pybind_pickle_class<mio::UncertainValue>(m, "UncertainValue")
        .def(py::init<ScalarType>(), py::arg("value") = 0.0)
        .def_property(
            "value",
            [](mio::UncertainValue& self) {
                return ScalarType(self);
            },
            [](mio::UncertainValue& self, ScalarType v) {
                self = v;
            })
        .def("set_distribution", //a property would be nicer but getter and setter use a different type
             &mio::UncertainValue::set_distribution)
        .def(
            "get_distribution",
            [](const mio::UncertainValue& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def(
            "get_distribution",
            [](mio::UncertainValue& self) {
                return self.get_distribution().get();
            },
            py::return_value_policy::reference_internal)
        .def("draw_sample", &mio::UncertainValue::draw_sample);

    auto contact_matrix_class = py::class_<mio::ContactMatrix>(m, "ContactMatrix");
    pymio::bind_damping_expression_members(contact_matrix_class);
    contact_matrix_class.def_property_readonly("num_groups", &mio::ContactMatrix::get_num_groups);

    auto contact_matrix_group_class = py::class_<mio::ContactMatrixGroup>(m, "ContactMatrixGroup");
    pymio::bind_damping_expression_group_members(contact_matrix_group_class);
    contact_matrix_group_class.def_property_readonly("num_groups", &mio::ContactMatrixGroup::get_num_groups);

    pymio::pybind_pickle_class<mio::DampingSampling>(m, "DampingSampling")
        .def(py::init([](const mio::UncertainValue& value, int level, int type, double time,
                         const std::vector<size_t>& matrices, const Eigen::Ref<const Eigen::VectorXd>& groups) {
                 return mio::DampingSampling(value, mio::DampingLevel(level), mio::DampingType(type),
                                             mio::SimulationTime(time), matrices, groups);
             }),
             py::arg("value"), py::arg("level"), py::arg("type"), py::arg("time"), py::arg("matrix_indices"),
             py::arg("group_weights"))
        .def_property("value", py::overload_cast<>(&mio::DampingSampling::get_value), &mio::DampingSampling::set_value,
                      py::return_value_policy::reference_internal)
        .def_property(
            "level",
            [](const mio::DampingSampling& self) {
                return int(self.get_level());
            },
            [](mio::DampingSampling& self, int lvl) {
                self.set_level(mio::DampingLevel(lvl));
            })
        .def_property(
            "type",
            [](const mio::DampingSampling& self) {
                return int(self.get_type());
            },
            [](mio::DampingSampling& self, int typ) {
                self.set_type(mio::DampingType(typ));
            })
        .def_property(
            "time",
            [](const mio::DampingSampling& self) {
                return double(self.get_time());
            },
            [](mio::DampingSampling& self, double t) {
                self.set_time(mio::SimulationTime(t));
            })
        .def_property("matrix_indices", &mio::DampingSampling::get_matrix_indices,
                      &mio::DampingSampling::set_matrix_indices)
        .def_property(
            "group_weights", &mio::DampingSampling::get_group_weights,
            [](mio::DampingSampling& self, const Eigen::Ref<const Eigen::VectorXd>& v) {
                self.set_group_weights(v);
            },
            py::return_value_policy::reference_internal);

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
