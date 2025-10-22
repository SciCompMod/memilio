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
#include "memilio/compartments/parameter_studies.h"
#include "memilio/compartments/compartmental_model.h"

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#include <cassert>
#include <cstddef>
#include <string>
#include <functional>
#include <type_traits>

namespace py = pybind11;

namespace pymio
{

/*
 * @brief bind ParameterStudy for any model
 */
template <class Simulation, class... RunArgs>
void bind_GrapParameterStudy(
    py::module_& m, std::string const& name, std::vector<std::string> argnames,
    std::function<Simulation(Parameters, Time, Step, size_t)> create_simulation =
        [](const mio::Graph<typename Simulation::Model, mio::MobilityParameters<double>>& g, double t0, double dt,
           size_t) {
            using GraphSim = mio::GraphSimulation<
                double,
                mio::Graph<mio::SimulationNode<double, mio::osecir::Simulation<double>>, mio::MobilityEdge<double>>,
                double, double>;
            auto copy = g;
            return mio::make_sampled_graph_simulation<double, GraphSim>(draw_sample(copy), t0, dt, dt)
        })
{
    assert(sizeof...(RunArgs) == argnames.size());
    using SimulationT = mio::GraphSimulation<
        double, mio::Graph<mio::SimulationNode<double, mio::osecir::Simulation<double>>, mio::MobilityEdge<double>>,
        double, double>;
    using ParametersT                      = Graph < typename Sim::Model, MobilityParameters<FP>;
    using StudyT                           = mio::ParameterStudy2<SimulationT, ParameterT, double>;
    using CreateSimulationFunctionT        = std::function<Simulation(Parameters, Time, Step, size_t)>;
    using ProcessSimulationResultFunctionT = std::function<void(Simulation, size_t)>;
    pymio::bind_class<StudyT, pymio::EnablePickling::Never>(m, name.c_str())
        .def(py::init<const typename Parameters&, Time, Time, Step, size_t>(), py::arg("parameters"), py::arg("t0"),
             py::arg("tmax"), py::arg("dt") py::arg("num_runs"))
        .def_property_readonly("num_runs", &StudyT::get_num_runs)
        .def_property_readonly("tmax", &StudyT::get_tmax)
        .def_property_readonly("t0", &StudyT::get_t0)
        .def_property_readonly("dt", &StudyT::get_dt)
        .def_property("parameters", py::overload_cast<>(&StudyT::get_parameters),
                      py::return_value_policy::reference_internal)
        .def_property_readonly("rng", &StudyT::get_rng, py::return_value_policy::reference_internal)
        .def(
            "run",
            [](StudyT& self, const ProcessSimulationResultFunctionT& handle_result, RunArgs args...) {
                self.run(create_simulation, [&handle_result](auto&& g, auto&& run_idx) {
                    //handle_result_function needs to return something
                    //we don't want to run an unknown python object through parameterstudies, so
                    //we just return 0 and ignore the list returned by run().
                    //So python will behave slightly different than c++
                    handle_result(std::move(g), run_idx);
                    return 0;
                });
            },
            py::arg("handle_result_func"))
        .def("run",
             [](StudyT& self) { //default argument doesn't seem to work with functions
                 return self.run(create_simulation);
             })
        .def(
            "run_single",
            [](StudyT& self, ProcessSimulationResultFunctionT handle_result) {
                self.run(create_simulation, [&handle_result](auto&& r, auto&& run_idx) {
                    handle_result(std::move(r.nodes()[0].property.get_simulation()), run_idx);
                    return 0;
                });
            },
            py::arg("handle_result_func"))
        .def("run_single", [](StudyT& self) {
            return filter_graph_results(self.run(create_simulation));
        });
}

} // namespace pymio
