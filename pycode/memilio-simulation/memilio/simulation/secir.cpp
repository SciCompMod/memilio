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
#include "secir/secir.h"
#include "secir/analyze_result.h"
#include "secir/parameter_space.h"
#include "memilio/compartments/parameter_studies.h"
#include "Eigen/Core"
#include "pybind11/stl_bind.h"
#include <vector>

namespace py = pybind11;

namespace
{

//select only the first node of the graph of each run, used for parameterstudy with single nodes
template <class Sim>
std::vector<Sim>
filter_graph_results(std::vector<mio::Graph<mio::SimulationNode<Sim>, mio::MigrationEdge>>&& graph_results)
{
    std::vector<Sim> results;
    results.reserve(graph_results.size());
    for (auto i = size_t(0); i < graph_results.size(); ++i) {
        results.emplace_back(std::move(graph_results[i].nodes()[0].property.get_simulation()));
    }
    return std::move(results);
}

/*
 * @brief bind ParameterStudy for any model
 */
template <class Simulation>
void bind_ParameterStudy(py::module& m, std::string const& name)
{
    py::class_<mio::ParameterStudy<Simulation>>(m, name.c_str())
        .def(py::init<const typename Simulation::Model&, double, double, size_t>(), py::arg("model"), py::arg("t0"),
             py::arg("tmax"), py::arg("num_runs"))
        .def(py::init<const mio::Graph<typename Simulation::Model, mio::MigrationParameters>&, double, double, double,
                      size_t>(),
             py::arg("model_graph"), py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("num_runs"))
        .def_property("num_runs", &mio::ParameterStudy<Simulation>::get_num_runs,
                      &mio::ParameterStudy<Simulation>::set_num_runs)
        .def_property("tmax", &mio::ParameterStudy<Simulation>::get_tmax, &mio::ParameterStudy<Simulation>::set_tmax)
        .def_property("t0", &mio::ParameterStudy<Simulation>::get_t0, &mio::ParameterStudy<Simulation>::set_t0)
        .def_property_readonly("model", py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_model),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("model", py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_model, py::const_),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_secir_model_graph),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_secir_model_graph, py::const_),
                               py::return_value_policy::reference_internal)
        .def(
            "run",
            [](mio::ParameterStudy<Simulation>& self,
               std::function<void(mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>)> handle_result) {
                self.run(
                    [](auto&& g) {
                        return draw_sample(g);
                    },
                    [&handle_result](auto&& g) {
                        handle_result(std::move(g));
                    });
            },
            py::arg("handle_result_func"))
        .def("run",
             [](mio::ParameterStudy<Simulation>& self) { //default argument doesn't seem to work with functions
                 return self.run([](auto&& g) {
                     return draw_sample(g);
                 });
             })
        .def(
            "run_single",
            [](mio::ParameterStudy<Simulation>& self, std::function<void(Simulation)> handle_result) {
                self.run(
                    [](auto&& g) {
                        return draw_sample(g);
                    },
                    [&handle_result](auto&& r) {
                        handle_result(std::move(r.nodes()[0].property.get_simulation()));
                    });
            },
            py::arg("handle_result_func"))
        .def("run_single", [](mio::ParameterStudy<Simulation>& self) {
            return filter_graph_results(self.run([](auto&& g) {
                return draw_sample(g);
            }));
        });
}

using Simulation = mio::SecirSimulation<>;
using MigrationGraph = mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>;

} // namespace

PYBIND11_MAKE_OPAQUE(std::vector<MigrationGraph>);

namespace pymio
{

//specialization of pretty_name for secir types
template <>
std::string pretty_name<mio::InfectionState>()
{
    return "InfectionState";
}
template <>
std::string pretty_name<mio::AgeGroup>()
{
    return "AgeGroup";
}

} // namespace pymio

PYBIND11_MODULE(_simulation_secir, m)
{
    // https://github.com/pybind/pybind11/issues/1153
    m.def("interpolate_simulation_result", static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&)>(
                                               &mio::interpolate_simulation_result));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<mio::TimeSeries<double>>);

    m.def("ensemble_mean", &mio::ensemble_mean);
    m.def("ensemble_percentile", &mio::ensemble_percentile);

    py::enum_<mio::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::InfectionState::Susceptible)
        .value("Exposed", mio::InfectionState::Exposed)
        .value("Carrier", mio::InfectionState::Carrier)
        .value("Infected", mio::InfectionState::Infected)
        .value("Hospitalized", mio::InfectionState::Hospitalized)
        .value("ICU", mio::InfectionState::ICU)
        .value("Recovered", mio::InfectionState::Recovered)
        .value("Dead", mio::InfectionState::Dead)
        .value("Count", mio::InfectionState::Count)
        .export_values();

    pymio::bind_Index<mio::InfectionState>(m, "Index_InfectionState");
    pymio::bind_Index<mio::AgeGroup>(m, "Index_AgeGroup");

    py::class_<mio::AgeGroup, mio::Index<mio::AgeGroup>>(m, "AgeGroup").def(py::init<size_t>());

    pymio::bind_MultiIndex<mio::AgeGroup, mio::InfectionState>(m, "Index_Agegroup_InfectionState");
    pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup, mio::InfectionState>(m, "SecirPopulationArray");
    pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, "AgeGroupArray");

    pymio::bind_ParameterSet<mio::SecirParamsBase>(m, "SecirParamsBase");

    py::class_<mio::SecirParams, mio::SecirParamsBase>(m, "SecirParams")
        .def(py::init<mio::AgeGroup>())
        .def("check_constraints", &mio::SecirParams::check_constraints)
        .def("apply_constraints", &mio::SecirParams::apply_constraints);

    pymio::bind_Population<mio::AgeGroup, mio::InfectionState>(m, "SecirPopulation");

    using SecirPopulations = mio::Populations<mio::AgeGroup, mio::InfectionState>;
    pymio::bind_CompartmentalModel<mio::InfectionState, SecirPopulations, mio::SecirParams>(m, "SecirModelBase");
    py::class_<mio::SecirModel, mio::CompartmentalModel<mio::InfectionState, SecirPopulations, mio::SecirParams>>(m, "SecirModel")
        .def(py::init<int>(), py::arg("num_agegroups"));

    pymio::bind_Simulation<mio::SecirSimulation<>>(m, "SecirSimulation");

    m.def("simulate", [](double t0, double tmax, double dt, const mio::SecirModel& model) { return mio::simulate(t0, tmax, dt, model); },
          "Simulates a SecirModel1 from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
          py::arg("model"));

    pymio::bind_SecirModelNode<mio::SecirModel>(m, "SecirModelNode");
    pymio::bind_SecirSimulationNode<mio::SecirSimulation<>>(m, "SecirSimulationNode");
    pymio::bind_SecirModelGraph<mio::SecirModel>(m, "SecirModelGraph");
    pymio::bind_MigrationGraph<Simulation>(m, "MigrationGraph");
    pymio::bind_GraphSimulation<MigrationGraph>(m, "MigrationSimulation");

    //normally, std::vector is bound to any python iterable, but this doesn't work for move-only elements
    //Bound the vector as a custom type that serves as output of ParameterStudy::run and input to
    //interpolate_ensemble_results
    py::bind_vector<std::vector<MigrationGraph>>(m, "EnsembleGraphResults");
    bind_ParameterStudy<mio::SecirSimulation<>>(m, "ParameterStudy");

    m.def("set_params_distributions_normal", &mio::set_params_distributions_normal, py::arg("model"), py::arg("t0"),
          py::arg("tmax"), py::arg("dev_rel"));

    m.def("draw_sample", [](mio::SecirModel& model) { return mio::draw_sample(model); }, py::arg("model"));

    m.def("interpolate_simulation_result",
          py::overload_cast<const MigrationGraph&>(&mio::interpolate_simulation_result<Simulation>));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<MigrationGraph>);

    m.attr("__version__") = "dev";
}
