/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef PYMIO_PARAMETER_STUDIES_H
#define PYMIO_PARAMETER_STUDIES_H

#include "memilio/compartments/parameter_studies.h"

namespace py = pybind11;

namespace pymio
{

/**
 * @brief Bind ParameterStudy for any parameter type. The run functions must be defined on the returned bind object.
 */
template <class Parameters, class Time, class Step>
auto bind_ParameterStudy(py::module_& m, std::string const& name)
{
    using StudyT = mio::ParameterStudy<Parameters, Time, Step>;

    auto study_object = bind_class<StudyT, pymio::EnablePickling::Never>(m, name.c_str());
    study_object
        .def(py::init<const Parameters&, Time, Time, Step, size_t>(), py::arg("parameters"), py::arg("t0"),
             py::arg("tmax"), py::arg("dt"), py::arg("num_runs"))
        .def_property_readonly("num_runs", &StudyT::get_num_runs)
        .def_property_readonly("tmax", &StudyT::get_tmax)
        .def_property_readonly("t0", &StudyT::get_t0)
        .def_property_readonly("dt", &StudyT::get_dt)
        .def_property_readonly("parameters", py::overload_cast<>(&StudyT::get_parameters),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("parameters", py::overload_cast<>(&StudyT::get_parameters, py::const_),
                               py::return_value_policy::reference_internal)
        .doc() = "Base binding of a ParameterStudy. Requires definitions for run method(s).";
    return study_object;
}

/**
 * @brief Bind ParameterStudy for ODE models (without using graphs).
 */
template <class Sim>
void bind_OdeParameterStudy(py::module_& m, std::string const& name)
{
    using ParametersT            = typename Sim::Model;
    using StudyT                 = mio::ParameterStudy<ParametersT, double>;
    const auto create_simulation = [](const ParametersT& m, double t0, double dt, size_t) {
        auto copy = m;
        draw_sample(copy);
        auto sim = Sim(std::move(copy), t0, dt);
        return sim;
    };

    bind_ParameterStudy<ParametersT, double, double>(m, name)
        .def(
            "run",
            [&create_simulation](StudyT& self, std::function<void(Sim, size_t)> handle_result) {
                self.run_serial(create_simulation, [&handle_result](auto&& g, auto&& run_idx) {
                    handle_result(std::move(g), run_idx);
                });
            },
            "Repeatedly simulate a model, handling its result in the provided function.", py::arg("handle_result_func"))
        .def(
            "run",
            [&create_simulation](StudyT& self) { //default argument doesn't seem to work with functions
                return self.run_serial(create_simulation, [](Sim&& result, size_t) {
                    return std::move(result);
                });
            },
            "Repeatedly simulate a model, returning all simulations in a vector.")
        .doc() = "Run a simulation multiple times, drawing new parameters each run.";
}

/**
 * @brief Bind ParameterStudy for ODE graph models.
 */
template <class Sim>
void bind_GraphOdeParameterStudy(py::module_& m, std::string const& name)
{
    using GraphT      = mio::Graph<mio::SimulationNode<double, Sim>, mio::MobilityEdge<double>>;
    using SimulationT = mio::GraphSimulation<double, GraphT, double, double>;
    using ParametersT = mio::Graph<typename Sim::Model, mio::MobilityParameters<double>>;
    using StudyT      = mio::ParameterStudy<ParametersT, double>;

    const auto create_simulation = [](const ParametersT& g, double t0, double dt, size_t) {
        auto copy = g;
        return mio::make_sampled_graph_simulation<double, Sim>(draw_sample(copy), t0, dt, dt);
    };

    bind_ParameterStudy<ParametersT, double, double>(m, name)
        .def(
            "run",
            [&create_simulation](StudyT& self, std::function<void(GraphT, size_t)> handle_result) {
                self.run_serial(create_simulation, [&handle_result](auto&& g, auto&& run_idx) {
                    handle_result(std::move(g.get_graph()), run_idx);
                });
            },
            "Repeatedly run a graph simulation, handling its result in the provided function.",
            py::arg("handle_result_func"))
        .def(
            "run",
            [&create_simulation](StudyT& self) { //default argument doesn't seem to work with functions
                return self.run_serial(create_simulation, [](SimulationT&& result, size_t) {
                    return std::move(result.get_graph());
                });
            },
            "Repeatedly run a graph simulation, returning all result graphs in a vector.")
        .doc() = "Run a graph simulation multiple times, drawing new parameters each run.";
}

} // namespace pymio

#endif // PYMIO_PARAMETER_STUDIES_H
