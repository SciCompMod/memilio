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
#include "utils/time_series.h"
#include "memilio/utils/time_series.h"
#include "pybind_util.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pybind11/eigen.h"

namespace py = pybind11;

namespace pymio
{

void bind_time_series(py::module_& m, std::string const& name)
{
    bind_class<mio::TimeSeries<double>, EnablePickling::Required>(m, name.c_str())
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
                    throw py::index_error("Index out of range."); //needs to throw exception for iterable
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
                    throw py::index_error("Index out of range."); //needs to throw exception for iterable
                }
            },
            py::is_operator(), py::arg("index"), py::arg("v"))
        .def(
            "print_table",
            [](const mio::TimeSeries<double>& self, bool return_string, const std::vector<std::string>& column_labels,
               size_t width, size_t precision, char separator, const std::string& header_prefix) {
                if (return_string) {
                    std::ostringstream oss;
                    self.print_table(oss, column_labels, width, precision, separator, header_prefix);
                    return py::object(py::str(oss.str()));
                }
                else {
                    self.print_table(column_labels, width, precision, separator, header_prefix);
                    return py::object(py::none());
                }
            },
            "Prints the TimeSeries as a formatted table. If return_string is true, the table is returned as a "
            "string. Otherwise, it is printed to the console.",
            py::arg("return_string") = false, py::arg("column_labels") = std::vector<std::string>{},
            py::arg("width") = 16, py::arg("precision") = 5, py::arg("separator") = ' ',
            py::arg("header_prefix") = "\n")
        .def(
            "export_csv",
            [](const mio::TimeSeries<double>& self, const std::string& filename,
               const std::vector<std::string>& column_labels, char separator, int precision) {
                auto result = self.export_csv(filename, column_labels, separator, precision);
                if (!result) {
                    throw py::value_error(result.error().message());
                }
            },
            py::arg("filename"), py::arg("column_labels") = std::vector<std::string>{}, py::arg("separator") = ',',
            py::arg("precision") = 6)
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
        .def(
            "as_ndarray",
            [](mio::TimeSeries<double>& self) {
                auto m = Eigen::Map<mio::TimeSeries<double>::Matrix>(self.data(), self.get_num_rows(),
                                                                     self.get_num_time_points());
                return Eigen::Ref<mio::TimeSeries<double>::Matrix>(m);
            },
            py::return_value_policy::reference_internal)

        .def("get_times", [](const mio::TimeSeries<double>& self) {
            auto times_range = self.get_times();
            std::vector<double> times(times_range.begin(), times_range.end());
            return times;
        });
}

} // namespace pymio
