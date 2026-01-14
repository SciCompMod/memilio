/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Daniel Abele, Maximilian Betz
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
#include "utils/parameter_set.h"
#include "compartments/simulation.h"
#include "compartments/flow_simulation.h"
#include "compartments/compartmental_model.h"
#include "compartments/parameter_studies.h"
#include "mobility/graph_simulation.h"
#include "mobility/metapopulation_mobility_instant.h"
#include "epidemiology/age_group.h"
#include "epidemiology/populations.h"
#include "io/mobility_io.h"
#include "io/result_io.h"

//Includes from MEmilio
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/analyze_result.h"
#include "ode_secirvvs/parameter_space.h"
#include "ode_secirvvs/parameters_io.h"
#include "memilio/data/analyze_result.h"
#include "memilio/compartments/flow_simulation.h"

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include <vector>

namespace py = pybind11;

namespace
{

enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

using MobilityGraph =
    mio::Graph<mio::SimulationNode<double, mio::osecirvvs::Simulation<double>>, mio::MobilityEdge<double>>;

} // namespace

namespace pymio
{
//specialization of pretty_name
template <>
inline std::string pretty_name<mio::osecirvvs::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MAKE_OPAQUE(std::vector<MobilityGraph>);

PYBIND11_MODULE(_simulation_osecirvvs, m)
{
    m.def("interpolate_simulation_result",
          static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&, const double)>(
              &mio::interpolate_simulation_result),
          py::arg("ts"), py::arg("abs_tol") = 1e-14);

    m.def("interpolate_simulation_result",
          static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&, const std::vector<double>&)>(
              &mio::interpolate_simulation_result),
          py::arg("ts"), py::arg("interpolation_times"));

    pymio::iterable_enum<mio::osecirvvs::InfectionState>(m, "InfectionState")
        .value("SusceptibleNaive", mio::osecirvvs::InfectionState::SusceptibleNaive)
        .value("SusceptiblePartialImmunity", mio::osecirvvs::InfectionState::SusceptiblePartialImmunity)
        .value("ExposedNaive", mio::osecirvvs::InfectionState::ExposedNaive)
        .value("ExposedPartialImmunity", mio::osecirvvs::InfectionState::ExposedPartialImmunity)
        .value("ExposedImprovedImmunity", mio::osecirvvs::InfectionState::ExposedImprovedImmunity)
        .value("InfectedNoSymptomsNaive", mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive)
        .value("InfectedNoSymptomsPartialImmunity", mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity)
        .value("InfectedNoSymptomsImprovedImmunity", mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity)
        .value("InfectedNoSymptomsNaiveConfirmed", mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed)
        .value("InfectedNoSymptomsPartialImmunityConfirmed",
               mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed)
        .value("InfectedNoSymptomsImprovedImmunityConfirmed",
               mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed)
        .value("InfectedSymptomsNaive", mio::osecirvvs::InfectionState::InfectedSymptomsNaive)
        .value("InfectedSymptomsPartialImmunity", mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity)
        .value("InfectedSymptomsImprovedImmunity", mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity)
        .value("InfectedSymptomsNaiveConfirmed", mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed)
        .value("InfectedSymptomsPartialImmunityConfirmed",
               mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed)
        .value("InfectedSymptomsImprovedImmunityConfirmed",
               mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed)
        .value("InfectedSevereNaive", mio::osecirvvs::InfectionState::InfectedSevereNaive)
        .value("InfectedSeverePartialImmunity", mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity)
        .value("InfectedSevereImprovedImmunity", mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity)
        .value("InfectedCriticalNaive", mio::osecirvvs::InfectionState::InfectedCriticalNaive)
        .value("InfectedCriticalPartialImmunity", mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity)
        .value("InfectedCriticalImprovedImmunity", mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity)
        .value("SusceptibleImprovedImmunity", mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity)
        .value("DeadNaive", mio::osecirvvs::InfectionState::DeadNaive)
        .value("DeadPartialImmunity", mio::osecirvvs::InfectionState::DeadPartialImmunity)
        .value("DeadImprovedImmunity", mio::osecirvvs::InfectionState::DeadImprovedImmunity);

    pymio::bind_ParameterSet<mio::osecirvvs::ParametersBase<double>, pymio::EnablePickling::Required>(m,
                                                                                                      "ParametersBase");

    pymio::bind_class<mio::osecirvvs::Parameters<double>, pymio::EnablePickling::Required,
                      mio::osecirvvs::ParametersBase<double>>(m, "Parameters")
        .def(py::init<mio::AgeGroup>())
        .def_property(
            "commuter_nondetection",
            [](const mio::osecirvvs::Parameters<double>& self) {
                return self.get_commuter_nondetection();
            },
            [](mio::osecirvvs::Parameters<double>& self, double v) {
                self.get_commuter_nondetection() = v;
            })
        .def_property(
            "start_commuter_detection",
            [](const mio::osecirvvs::Parameters<double>& self) {
                return self.get_start_commuter_detection();
            },
            [](mio::osecirvvs::Parameters<double>& self, double v) {
                self.get_start_commuter_detection() = v;
            })
        .def_property(
            "end_commuter_detection",
            [](const mio::osecirvvs::Parameters<double>& self) {
                return self.get_end_commuter_detection();
            },
            [](mio::osecirvvs::Parameters<double>& self, double v) {
                self.get_end_commuter_detection() = v;
            })
        .def_property(
            "end_dynamic_npis",
            [](const mio::osecirvvs::Parameters<double>& self) {
                return self.get_end_dynamic_npis();
            },
            [](mio::osecirvvs::Parameters<double>& self, double v) {
                self.get_end_dynamic_npis() = v;
            })
        .def("check_constraints", &mio::osecirvvs::Parameters<double>::check_constraints)
        .def("apply_constraints", &mio::osecirvvs::Parameters<double>::apply_constraints);

    using Populations = mio::Populations<double, mio::AgeGroup, mio::osecirvvs::InfectionState>;
    pymio::bind_Population(m, "Populations", mio::Tag<mio::osecirvvs::Model<double>::Populations>{});

    pymio::bind_CompartmentalModel<mio::osecirvvs::InfectionState, Populations, mio::osecirvvs::Parameters<double>,
                                   pymio::EnablePickling::Never>(m, "ModelBase");
    pymio::bind_class<mio::osecirvvs::Model<double>, pymio::EnablePickling::Required,
                      mio::CompartmentalModel<double, mio::osecirvvs::InfectionState, Populations,
                                              mio::osecirvvs::Parameters<double>>>(m, "Model")
        .def(py::init<int>(), py::arg("num_agegroups"));

    pymio::bind_Simulation<mio::osecirvvs::Simulation<double>>(m, "Simulation");
    pymio::bind_Flow_Simulation<
        mio::osecirvvs::Simulation<double, mio::FlowSimulation<double, mio::osecirvvs::Model<double>>>>(
        m, "FlowSimulation");

    m.def("simulate", &mio::osecirvvs::simulate<double>, "Simulates an ODE SECIRVVS model from t0 to tmax.",
          py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none());

    m.def("simulate_flows", &mio::osecirvvs::simulate_flows<double>,
          "Simulates an ODE SECIRVVS model with flows from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
          py::arg("model"), py::arg("integrator") = py::none());

    pymio::bind_ModelNode<mio::osecirvvs::Model<double>>(m, "ModelNode");
    pymio::bind_SimulationNode<mio::osecirvvs::Simulation<double>>(m, "SimulationNode");
    pymio::bind_ModelGraph<mio::osecirvvs::Model<double>>(m, "ModelGraph");
    pymio::bind_MobilityGraph<mio::osecirvvs::Simulation<double>>(m, "MobilityGraph");
    pymio::bind_GraphSimulation<MobilityGraph>(m, "MobilitySimulation");

    //normally, std::vector is bound to any python iterable, but this doesn't work for move-only elements
    //Bound the vector as a custom type that serves as output of ParameterStudy::run and input to
    //interpolate_ensemble_results
    py::bind_vector<std::vector<MobilityGraph>>(m, "EnsembleGraphResults");
    pymio::bind_OdeParameterStudy<mio::osecirvvs::Simulation<double>>(m, "ParameterStudy");
    pymio::bind_GraphOdeParameterStudy<mio::osecirvvs::Simulation<double>>(m, "GraphParameterStudy");

    m.def(
        "draw_sample",
        [](mio::osecirvvs::Model<double>& model) {
            return mio::osecirvvs::draw_sample(model);
        },
        py::arg("model"));

    // These functions are in general not secirvvs dependent, only with the current config
    m.def(
        "set_nodes",
        [](const mio::osecirvvs::Parameters<double>& params, mio::Date start_date, mio::Date end_date,
           const std::string& data_dir, const std::string& population_data_path, bool is_node_for_county,
           mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>>& params_graph,
           const std::vector<double>& scaling_factor_inf, double scaling_factor_icu, double tnt_capacity_factor,
           int num_days = 0, bool export_time_series = false) {
            auto result =
                mio::set_nodes<double, // FP
                               mio::osecirvvs::TestAndTraceCapacity<double>, mio::osecirvvs::ContactPatterns<double>,
                               mio::osecirvvs::Model<double>, mio::MobilityParameters<double>,
                               mio::osecirvvs::Parameters<double>,
                               decltype(mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model<ScalarType>>),
                               decltype(mio::get_node_ids)>(
                    params, start_date, end_date, data_dir, population_data_path, is_node_for_county, params_graph,
                    mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model<double>>, mio::get_node_ids,
                    scaling_factor_inf, scaling_factor_icu, tnt_capacity_factor, num_days, export_time_series);
            return pymio::check_and_throw(result);
        },
        py::return_value_policy::move);

    m.def(
        "set_edges",
        [](const std::string& mobility_data_file,
           mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>>& params_graph,
           size_t contact_locations_size) {
            auto mobile_comp = {mio::osecirvvs::InfectionState::SusceptibleNaive,
                                mio::osecirvvs::InfectionState::ExposedNaive,
                                mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive,
                                mio::osecirvvs::InfectionState::InfectedSymptomsNaive,
                                mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity,
                                mio::osecirvvs::InfectionState::SusceptiblePartialImmunity,
                                mio::osecirvvs::InfectionState::ExposedPartialImmunity,
                                mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity,
                                mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity,
                                mio::osecirvvs::InfectionState::ExposedImprovedImmunity,
                                mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity,
                                mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity};
            auto weights     = std::vector<ScalarType>{0., 0., 1.0, 1.0, 0.33, 0., 0.};
            auto result      = mio::set_edges<double, // FP,
                                              ContactLocation, mio::osecirvvs::Model<double>,
                                              mio::MobilityParameters<double>, mio::MobilityCoefficientGroup<double>,
                                              mio::osecirvvs::InfectionState, decltype(mio::read_mobility_plain)>(
                mobility_data_file, params_graph, mobile_comp, contact_locations_size, mio::read_mobility_plain,
                weights);
            return pymio::check_and_throw(result);
        },
        py::return_value_policy::move);

#ifdef MEMILIO_HAS_HDF5
    pymio::bind_save_results<mio::osecirvvs::Model<double>>(m);
#endif // MEMILIO_HAS_HDF5

#ifdef MEMILIO_HAS_JSONCPP
    pymio::bind_write_graph<mio::osecirvvs::Model<double>>(m);
    pymio::bind_read_graph<mio::osecirvvs::Model<double>>(m);
    m.def(
        "read_input_data_county",
        [](std::vector<mio::osecirvvs::Model<double>>& model, mio::Date date, const std::vector<int>& county,
           const std::vector<double>& scaling_factor_inf, double scaling_factor_icu, const std::string& dir,
           int num_days = 0, bool export_time_series = false) {
            auto result = mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model<double>>(
                model, date, county, scaling_factor_inf, scaling_factor_icu, dir, num_days, export_time_series);
            return pymio::check_and_throw(result);
        },
        py::return_value_policy::move);
#endif // MEMILIO_HAS_JSONCPP

    m.def("interpolate_simulation_result",
          py::overload_cast<const MobilityGraph&>(
              &mio::interpolate_simulation_result<double, mio::osecirvvs::Simulation<double>>));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<MobilityGraph>);

    m.attr("__version__") = "dev";
}
