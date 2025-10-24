/* 
* Copyright (C) 2020-2025 MEmilio
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

//Includes from pymio
#include "pybind_util.h"
#include "compartments/simulation.h"
#include "compartments/flow_simulation.h"
#include "compartments/compartmental_model.h"
#include "epidemiology/age_group.h"
#include "epidemiology/populations.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "utils/index.h"
#include "mobility/graph_simulation.h"
#include "mobility/metapopulation_mobility_instant.h"
#include "io/mobility_io.h"
#include "io/result_io.h"

//Includes from MEmilio
#include "ode_secir/model.h"
#include "ode_secir/analyze_result.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/parameters_io.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/data/analyze_result.h"
#include "memilio/mobility/graph.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/epi_data.h"
#include "memilio/config.h"

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "Eigen/Core"
#include <algorithm>
#include <cstddef>
#include <vector>

namespace py = pybind11;

namespace
{

/*
 * @brief bind ParameterStudy for any model
 */
template <class Sim>
void bind_ParameterStudy(py::module_& m, std::string const& name)
{
    using GraphT      = mio::Graph<mio::SimulationNode<double, Sim>, mio::MobilityEdge<double>>;
    using SimulationT = mio::GraphSimulation<double, GraphT, double, double>;
    using ParametersT = mio::Graph<typename Sim::Model, mio::MobilityParameters<double>>;
    using StudyT      = mio::ParameterStudy2<ParametersT, double>;

    const auto create_simulation = [](const ParametersT& g, double t0, double dt, size_t) {
        auto copy = g;
        return mio::make_sampled_graph_simulation<double, Sim>(draw_sample(copy), t0, dt, dt);
    };

    pymio::bind_class<StudyT, pymio::EnablePickling::Never>(m, name.c_str())
        .def(py::init<const ParametersT&, double, double, double, size_t>(), py::arg("parameters"), py::arg("t0"),
             py::arg("tmax"), py::arg("dt"), py::arg("num_runs"))
        .def_property_readonly("num_runs", &StudyT::get_num_runs)
        .def_property_readonly("tmax", &StudyT::get_tmax)
        .def_property_readonly("t0", &StudyT::get_t0)
        .def_property_readonly("dt", &StudyT::get_dt)
        .def_property_readonly("model_graph", py::overload_cast<>(&StudyT::get_parameters),
                               py::return_value_policy::reference_internal)
        .def_property_readonly("model_graph", py::overload_cast<>(&StudyT::get_parameters, py::const_),
                               py::return_value_policy::reference_internal)
        .def(
            "run",
            [&create_simulation](StudyT& self, std::function<void(GraphT, size_t)> handle_result) {
                self.run_serial(create_simulation, [&handle_result](auto&& g, auto&& run_idx) {
                    //handle_result_function needs to return something
                    //we don't want to run an unknown python object through parameterstudies, so
                    //we just return 0 and ignore the list returned by run().
                    //So python will behave slightly different than c++
                    handle_result(std::move(g.get_graph()), run_idx);
                    return 0;
                });
            },
            py::arg("handle_result_func"))
        .def("run",
             [&create_simulation](StudyT& self) { //default argument doesn't seem to work with functions
                 return self.run_serial(create_simulation, [](SimulationT&& result, size_t) {
                     return std::move(result.get_graph());
                 });
             })
        .def(
            "run_single",
            [&create_simulation](StudyT& self, std::function<void(Sim, size_t)> handle_result) {
                self.run_serial(create_simulation, [&handle_result](auto&& r, auto&& run_idx) {
                    handle_result(std::move(r.get_graph().nodes()[0].property.get_simulation()), run_idx);
                    return 0;
                });
            },
            py::arg("handle_result_func"))
        .def("run_single", [&create_simulation](StudyT& self) {
            return self.run_serial(create_simulation, [](SimulationT&& result, size_t) {
                //select only the first node of the graph of each run, used for parameterstudy with single nodes
                return std::move(result.get_graph().nodes()[0].property.get_simulation());
            });
        });
}

enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

using MobilityGraph =
    mio::Graph<mio::SimulationNode<double, mio::osecir::Simulation<double>>, mio::MobilityEdge<double>>;

} // namespace

namespace pymio
{
//specialization of pretty_name
template <>
inline std::string pretty_name<mio::osecir::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MAKE_OPAQUE(std::vector<MobilityGraph>);

PYBIND11_MODULE(_simulation_osecir, m)
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

    m.def("ensemble_mean", &mio::ensemble_mean<double>);
    m.def("ensemble_percentile", &mio::ensemble_percentile<double>);

    pymio::iterable_enum<mio::osecir::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::osecir::InfectionState::Susceptible)
        .value("Exposed", mio::osecir::InfectionState::Exposed)
        .value("InfectedNoSymptoms", mio::osecir::InfectionState::InfectedNoSymptoms)
        .value("InfectedNoSymptomsConfirmed", mio::osecir::InfectionState::InfectedNoSymptomsConfirmed)
        .value("InfectedSymptoms", mio::osecir::InfectionState::InfectedSymptoms)
        .value("InfectedSymptomsConfirmed", mio::osecir::InfectionState::InfectedSymptomsConfirmed)
        .value("InfectedSevere", mio::osecir::InfectionState::InfectedSevere)
        .value("InfectedCritical", mio::osecir::InfectionState::InfectedCritical)
        .value("Recovered", mio::osecir::InfectionState::Recovered)
        .value("Dead", mio::osecir::InfectionState::Dead);

    pymio::bind_ParameterSet<mio::osecir::ParametersBase<double>, pymio::EnablePickling::Required>(m, "ParametersBase");

    pymio::bind_class<mio::osecir::Parameters<double>, pymio::EnablePickling::Required,
                      mio::osecir::ParametersBase<double>>(m, "Parameters")
        .def(py::init<mio::AgeGroup>())
        .def("check_constraints", &mio::osecir::Parameters<double>::check_constraints)
        .def("apply_constraints", &mio::osecir::Parameters<double>::apply_constraints);

    using Populations = mio::Populations<double, mio::AgeGroup, mio::osecir::InfectionState>;
    pymio::bind_Population(m, "Populations", mio::Tag<mio::osecir::Model<double>::Populations>{});
    pymio::bind_CompartmentalModel<mio::osecir::InfectionState, Populations, mio::osecir::Parameters<double>,
                                   pymio::EnablePickling::Never>(m, "ModelBase");
    pymio::bind_class<
        mio::osecir::Model<double>, pymio::EnablePickling::Required,
        mio::CompartmentalModel<double, mio::osecir::InfectionState, Populations, mio::osecir::Parameters<double>>>(
        m, "Model")
        .def(py::init<int>(), py::arg("num_agegroups"));

    pymio::bind_Simulation<mio::osecir::Simulation<double>>(m, "Simulation");
    pymio::bind_Flow_Simulation<
        mio::osecir::Simulation<double, mio::FlowSimulation<double, mio::osecir::Model<double>>>>(m, "FlowSimulation");

    m.def("simulate", &mio::osecir::simulate<double>, "Simulates an ODE SECIHURD model from t0 to tmax.", py::arg("t0"),
          py::arg("tmax"), py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none());

    m.def("simulate_flows", &mio::osecir::simulate_flows<double>,
          "Simulates an ODE SECIHURD model with flows from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
          py::arg("model"), py::arg("integrator") = py::none());

    pymio::bind_ModelNode<mio::osecir::Model<double>>(m, "ModelNode");
    pymio::bind_SimulationNode<mio::osecir::Simulation<double>>(m, "SimulationNode");
    pymio::bind_ModelGraph<mio::osecir::Model<double>>(m, "ModelGraph");
    pymio::bind_MobilityGraph<mio::osecir::Simulation<double>>(m, "MobilityGraph");
    pymio::bind_GraphSimulation<MobilityGraph>(m, "MobilitySimulation");

    //normally, std::vector is bound to any python iterable, but this doesn't work for move-only elements
    //Bound the vector as a custom type that serves as output of ParameterStudy::run and input to
    //interpolate_ensemble_results
    py::bind_vector<std::vector<MobilityGraph>>(m, "EnsembleGraphResults");
    bind_ParameterStudy<mio::osecir::Simulation<double>>(m, "ParameterStudy");

    m.def("set_params_distributions_normal", &mio::osecir::set_params_distributions_normal<double>, py::arg("model"),
          py::arg("t0"), py::arg("tmax"), py::arg("dev_rel"));

    m.def(
        "draw_sample",
        [](mio::osecir::Model<double>& model) {
            return mio::osecir::draw_sample(model);
        },
        py::arg("model"));

    // These functions are in general not model dependent, only with the current config
    m.def(
        "set_nodes",
        [](const mio::osecir::Parameters<double>& params, mio::Date start_date, mio::Date end_date,
           const std::string& data_dir, const std::string& population_data_path, bool is_node_for_county,
           mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>& params_graph,
           const std::vector<double>& scaling_factor_inf, double scaling_factor_icu, double tnt_capacity_factor,
           int num_days = 0, bool export_time_series = false) {
            auto result = mio::set_nodes<
                double, // FP
                mio::osecir::TestAndTraceCapacity<double>, mio::osecir::ContactPatterns<double>,
                mio::osecir::Model<double>, mio::MobilityParameters<double>, mio::osecir::Parameters<double>,
                decltype(mio::osecir::read_input_data_county<mio::osecir::Model<double>>), decltype(mio::get_node_ids)>(
                params, start_date, end_date, data_dir, population_data_path, is_node_for_county, params_graph,
                mio::osecir::read_input_data_county<mio::osecir::Model<double>>, mio::get_node_ids, scaling_factor_inf,
                scaling_factor_icu, tnt_capacity_factor, num_days, export_time_series);
            return pymio::check_and_throw(result);
        },
        py::return_value_policy::move);

    pymio::iterable_enum<ContactLocation>(m, "ContactLocation")
        .value("Home", ContactLocation::Home)
        .value("School", ContactLocation::School)
        .value("Work", ContactLocation::Work)
        .value("Other", ContactLocation::Other);

    m.def(
        "set_edges",
        [](const std::string& mobility_data_file,
           mio::Graph<mio::osecir::Model<double>, mio::MobilityParameters<double>>& params_graph,
           size_t contact_locations_size) {
            auto mobile_comp = {mio::osecir::InfectionState::Susceptible, mio::osecir::InfectionState::Exposed,
                                mio::osecir::InfectionState::InfectedNoSymptoms,
                                mio::osecir::InfectionState::InfectedSymptoms, mio::osecir::InfectionState::Recovered};
            auto weights     = std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.};
            auto result      = mio::set_edges<double, // FP
                                              ContactLocation, mio::osecir::Model<double>, mio::MobilityParameters<double>,
                                              mio::MobilityCoefficientGroup<double>, mio::osecir::InfectionState,
                                              decltype(mio::read_mobility_plain)>(mobility_data_file, params_graph,
                                                                                  mobile_comp, contact_locations_size,
                                                                                  mio::read_mobility_plain, weights);
            return pymio::check_and_throw(result);
        },
        py::return_value_policy::move);

#ifdef MEMILIO_HAS_HDF5
    pymio::bind_save_results<mio::osecir::Model<double>>(m);
#endif // MEMILIO_HAS_HDF5

#ifdef MEMILIO_HAS_JSONCPP
    pymio::bind_write_graph<mio::osecir::Model<double>>(m);
    pymio::bind_read_graph<mio::osecir::Model<double>>(m);
    m.def(
        "read_input_data_county",
        [](std::vector<mio::osecir::Model<double>>& model, mio::Date date, const std::vector<int>& county,
           const std::vector<double>& scaling_factor_inf, double scaling_factor_icu, const std::string& dir,
           int num_days = 0, bool export_time_series = false) {
            auto result = mio::osecir::read_input_data_county<mio::osecir::Model<double>>(
                model, date, county, scaling_factor_inf, scaling_factor_icu, dir, num_days, export_time_series);
            return pymio::check_and_throw(result);
        },
        py::return_value_policy::move);
#endif // MEMILIO_HAS_JSONCPP

    m.def("interpolate_simulation_result",
          py::overload_cast<const MobilityGraph&>(
              &mio::interpolate_simulation_result<double, mio::osecir::Simulation<double>>));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<MobilityGraph>);

    m.attr("__version__") = "dev";
}
