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
#include "compartments/simulation.h"
#include "compartments/compartmentalmodel.h"
#include "epidemiology/populations.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "utils/index.h"
#include "mobility/graph_simulation.h"
#include "mobility/mobility.h"
#include "secir/secir.h"
#include "secir/parameter_studies.h"
#include "secir/analyze_result.h"
#include "Eigen/Core"
#include "pybind11/stl_bind.h"
#include <vector>

using namespace pybind11;

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
void bind_ParameterStudy(pybind11::module& m, std::string const& name)
{
    pybind11::class_<mio::ParameterStudy<Simulation>>(m, name.c_str())
        .def(pybind11::init<const typename Simulation::Model&, double, double, size_t>(), pybind11::arg("model"), pybind11::arg("t0"),
             pybind11::arg("tmax"), pybind11::arg("num_runs"))
        .def(pybind11::init<const typename Simulation::Model&, double, double, double, size_t>(), pybind11::arg("model"),
             pybind11::arg("t0"), pybind11::arg("tmax"), pybind11::arg("dev_rel"), pybind11::arg("num_runs"))
        .def(pybind11::init<const mio::Graph<typename Simulation::Model, mio::MigrationParameters>&, double, double, double,
                      size_t>(),
             pybind11::arg("model_graph"), pybind11::arg("t0"), pybind11::arg("tmax"), pybind11::arg("dt"), pybind11::arg("num_runs"))
        .def_property("num_runs", &mio::ParameterStudy<Simulation>::get_num_runs,
                      &mio::ParameterStudy<Simulation>::set_num_runs)
        .def_property("tmax", &mio::ParameterStudy<Simulation>::get_tmax, &mio::ParameterStudy<Simulation>::set_tmax)
        .def_property("t0", &mio::ParameterStudy<Simulation>::get_t0, &mio::ParameterStudy<Simulation>::set_t0)
        .def_property_readonly("model", pybind11::overload_cast<>(&mio::ParameterStudy<Simulation>::get_model),
                               pybind11::return_value_policy::reference_internal)
        .def_property_readonly("model", pybind11::overload_cast<>(&mio::ParameterStudy<Simulation>::get_model, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               pybind11::overload_cast<>(&mio::ParameterStudy<Simulation>::get_secir_model_graph),
                               pybind11::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               pybind11::overload_cast<>(&mio::ParameterStudy<Simulation>::get_secir_model_graph, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def(
            "run",
            [](mio::ParameterStudy<Simulation>& self, std::function<void(mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>)> handle_result) {
                self.run([&handle_result](auto&& g) { handle_result(std::move(g)); });
            },
            pybind11::arg("handle_result_func"))
        .def("run",
             [](mio::ParameterStudy<Simulation>& self) { //default argument doesn't seem to work with functions
                 return self.run();
             })
        .def(
            "run_single",
            [](mio::ParameterStudy<Simulation>& self, std::function<void(Simulation)> handle_result) {
                self.run([&handle_result](auto&& r) {
                    handle_result(std::move(r.nodes()[0].property.get_simulation()));
                });
            },
            pybind11::arg("handle_result_func"))
        .def("run_single", [](mio::ParameterStudy<Simulation>& self) {
            return filter_graph_results(self.run());
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

    pymio::iterable_enum<mio::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::InfectionState::Susceptible)
        .value("Exposed", mio::InfectionState::Exposed)
        .value("Carrier", mio::InfectionState::Carrier)
        .value("Infected", mio::InfectionState::Infected)
        .value("Hospitalized", mio::InfectionState::Hospitalized)
        .value("ICU", mio::InfectionState::ICU)
        .value("Recovered", mio::InfectionState::Recovered)
        .value("Dead", mio::InfectionState::Dead);

    pymio::bind_Index<mio::InfectionState>(m, "Index_InfectionState");
    pymio::bind_Index<mio::AgeGroup>(m, "Index_AgeGroup");

    pybind11::class_<mio::AgeGroup, mio::Index<mio::AgeGroup>>(m, "AgeGroup").def(pybind11::init<size_t>());

    pymio::bind_MultiIndex<mio::AgeGroup, mio::InfectionState>(m, "Index_Agegroup_InfectionState");
    pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup, mio::InfectionState>(m, "SecirPopulationArray");
    pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, "AgeGroupArray");

    pymio::bind_ParameterSet<mio::SecirParamsBase>(m, "SecirParamsBase");

    pybind11::class_<mio::SecirParams, mio::SecirParamsBase>(m, "SecirParams")
        .def(pybind11::init<mio::AgeGroup>())
        .def("check_constraints", &mio::SecirParams::check_constraints)
        .def("apply_constraints", &mio::SecirParams::apply_constraints);

    pymio::bind_Population<mio::AgeGroup, mio::InfectionState>(m, "SecirPopulation");

    using SecirPopulations = mio::Populations<mio::AgeGroup, mio::InfectionState>;
    pymio::bind_CompartmentalModel<SecirPopulations, mio::SecirParams>(m, "SecirModelBase");
    pybind11::class_<mio::SecirModel, mio::CompartmentalModel<SecirPopulations, mio::SecirParams>>(m, "SecirModel")
        .def(pybind11::init<int>(), pybind11::arg("num_agegroups"));

    pymio::bind_Simulation<mio::SecirSimulation<>>(m, "SecirSimulation");

    m.def("simulate", [](double t0, double tmax, double dt, const mio::SecirModel& model) { return mio::simulate(t0, tmax, dt, model); },
          "Simulates a SecirModel1 from t0 to tmax.", pybind11::arg("t0"), pybind11::arg("tmax"), pybind11::arg("dt"),
          pybind11::arg("model"));

    pymio::bind_ModelNode<mio::SecirModel>(m, "SecirModelNode");
    pymio::bind_SimulationNode<mio::SecirSimulation<>>(m, "SecirSimulationNode");
    pymio::bind_ModelGraph<mio::SecirModel>(m, "SecirModelGraph");
    pymio::bind_MigrationGraph<Simulation>(m, "MigrationGraph");
    pymio::bind_GraphSimulation<MigrationGraph>(m, "MigrationSimulation");

    //normally, std::vector is bound to any python iterable, but this doesn't work for move-only elements
    //Bound the vector as a custom type that serves as output of ParameterStudy::run and input to
    //interpolate_ensemble_results
    pybind11::bind_vector<std::vector<MigrationGraph>>(m, "EnsembleGraphResults");
    bind_ParameterStudy<mio::SecirSimulation<>>(m, "ParameterStudy");

    m.def("set_params_distributions_normal", &mio::set_params_distributions_normal, pybind11::arg("model"), pybind11::arg("t0"),
          pybind11::arg("tmax"), pybind11::arg("dev_rel"));

    m.def("draw_sample", &mio::draw_sample, pybind11::arg("model"));

    m.def("interpolate_simulation_result",
          pybind11::overload_cast<const MigrationGraph&>(&mio::interpolate_simulation_result<Simulation>));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<MigrationGraph>);

    m.attr("__version__") = "dev";
}
