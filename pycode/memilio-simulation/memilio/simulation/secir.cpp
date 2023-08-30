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
#include "mobility/metapopulation_mobility_instant.h"
#include "ode_secir/model.h"
#include "ode_secir/analyze_result.h"
#include "ode_secir/parameter_space.h"
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
void bind_ParameterStudy(py::module_& m, std::string const& name)
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
                               py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_model_graph),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("secir_model_graph",
                               py::overload_cast<>(&mio::ParameterStudy<Simulation>::get_model_graph, py::const_),
                               py::return_value_policy::reference_internal)
        .def("run",
            [](mio::ParameterStudy<Simulation>& self,
               std::function<void(mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>, size_t)> handle_result) {
                self.run(
                    [](auto&& g) {
                        return draw_sample(g);
                    },
                    [&handle_result](auto&& g, auto&& run_idx) {
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
             [](mio::ParameterStudy<Simulation>& self) { //default argument doesn't seem to work with functions
                 return self.run([](auto&& g) {
                     return draw_sample(g);
                 });
             })
        .def(
            "run_single",
            [](mio::ParameterStudy<Simulation>& self, std::function<void(Simulation, size_t)> handle_result) {
                self.run(
                    [](auto&& g) {
                        return draw_sample(g);
                    },
                    [&handle_result](auto&& r, auto&& run_idx) {
                        handle_result(std::move(r.nodes()[0].property.get_simulation()), run_idx);
                        return 0;
                    });
            },
            py::arg("handle_result_func"))
        .def("run_single", [](mio::ParameterStudy<Simulation>& self) {
            return filter_graph_results(self.run([](auto&& g) {
                return draw_sample(g);
            }));
        });
}

using Simulation     = mio::osecir::Simulation<>;
using MigrationGraph = mio::Graph<mio::SimulationNode<Simulation>, mio::MigrationEdge>;

} // namespace

PYBIND11_MAKE_OPAQUE(std::vector<MigrationGraph>);

namespace pymio
{

//specialization of pretty_name for secir types
template <>
std::string pretty_name<mio::osecir::InfectionState>()
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
    m.def("interpolate_simulation_result",
          static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&, const double)>(
              &mio::interpolate_simulation_result),
          py::arg("ts"), py::arg("abs_tol") = 1e-14);

    m.def("interpolate_simulation_result",
          static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&, const std::vector<double>&)>(
              &mio::interpolate_simulation_result),
          py::arg("ts"), py::arg("interpolation_times"));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<mio::TimeSeries<double>>);

    m.def("ensemble_mean", &mio::ensemble_mean);
    m.def("ensemble_percentile", &mio::ensemble_percentile);

    pymio::iterable_enum<mio::osecir::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::osecir::InfectionState::Susceptible)
        .value("Exposed", mio::osecir::InfectionState::Exposed)
        .value("InfectedNoSymptoms", mio::osecir::InfectionState::InfectedNoSymptoms)
        .value("InfectedSymptoms", mio::osecir::InfectionState::InfectedSymptoms)
        .value("InfectedSevere", mio::osecir::InfectionState::InfectedSevere)
        .value("InfectedCritical", mio::osecir::InfectionState::InfectedCritical)
        .value("Recovered", mio::osecir::InfectionState::Recovered)
        .value("Dead", mio::osecir::InfectionState::Dead);

    pymio::bind_CustomIndexArray<mio::UncertainValue, mio::AgeGroup>(m, "AgeGroupArray");

    pymio::bind_ParameterSet<mio::osecir::ParametersBase>(m, "ParametersBase");

    py::class_<mio::osecir::Parameters, mio::osecir::ParametersBase>(m, "Parameters")
        .def(py::init<mio::AgeGroup>())
        .def("check_constraints", &mio::osecir::Parameters::check_constraints)
        .def("apply_constraints", &mio::osecir::Parameters::apply_constraints);

    using SecirPopulations = mio::Populations<mio::AgeGroup, mio::osecir::InfectionState>;
    pymio::bind_Population(m, "SecirPopulation", mio::Tag<mio::osecir::Model::Populations>{});
    py::class_<mio::AgeGroup, mio::Index<mio::AgeGroup>>(m, "AgeGroup").def(py::init<size_t>());
    pymio::bind_CompartmentalModel<mio::osecir::InfectionState, SecirPopulations, mio::osecir::Parameters>(m,
                                                                                                           "ModelBase");
    py::class_<mio::osecir::Model,
               mio::CompartmentalModel<mio::osecir::InfectionState, SecirPopulations, mio::osecir::Parameters>>(m,
                                                                                                                "Model")
        .def(py::init<int>(), py::arg("num_agegroups"));

    pymio::bind_Simulation<mio::osecir::Simulation<>>(m, "Simulation");

    m.def(
        "simulate",
        [](double t0, double tmax, double dt, const mio::osecir::Model& model) {
            return mio::simulate(t0, tmax, dt, model);
        },
        "Simulates a Secir Model1 from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"));

    pymio::bind_ModelNode<mio::osecir::Model>(m, "ModelNode");
    pymio::bind_SimulationNode<mio::osecir::Simulation<>>(m, "SimulationNode");
    pymio::bind_ModelGraph<mio::osecir::Model>(m, "ModelGraph");
    pymio::bind_MigrationGraph<Simulation>(m, "MigrationGraph");
    pymio::bind_GraphSimulation<MigrationGraph>(m, "MigrationSimulation");

    //normally, std::vector is bound to any python iterable, but this doesn't work for move-only elements
    //Bound the vector as a custom type that serves as output of ParameterStudy::run and input to
    //interpolate_ensemble_results
    py::bind_vector<std::vector<MigrationGraph>>(m, "EnsembleGraphResults");
    bind_ParameterStudy<mio::osecir::Simulation<>>(m, "ParameterStudy");

    m.def("set_params_distributions_normal", &mio::osecir::set_params_distributions_normal, py::arg("model"),
          py::arg("t0"), py::arg("tmax"), py::arg("dev_rel"));

    m.def(
        "draw_sample",
        [](mio::osecir::Model& model) {
            return mio::osecir::draw_sample(model);
        },
        py::arg("model"));

    m.def("interpolate_simulation_result",
          py::overload_cast<const MigrationGraph&>(&mio::interpolate_simulation_result<Simulation>));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<MigrationGraph>);

    m.attr("__version__") = "dev";
}
