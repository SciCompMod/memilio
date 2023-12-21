/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "pybind_util.h"

//Includes from pymio
#include "utils/parameter_set.h"
#include "compartments/simulation.h"
#include "compartments/compartmentalmodel.h"
#include "mobility/graph_simulation.h"
#include "mobility/metapopulation_mobility_instant.h"
#include "epidemiology/populations.h"
#include "io/mobility_io.h"

//Includes from MEmilio
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/analyze_result.h"
#include "ode_secirvvs/parameter_space.h"
#include "memilio/compartments/flow_simulation.h"

#include "pybind11/stl_bind.h"
#include <vector>

namespace py = pybind11;

namespace pymio
{
//specialization of pretty_name
template <>
std::string pretty_name<mio::AgeGroup>()
{
    return "AgeGroup";
}

template <>
std::string pretty_name<mio::osecirvvs::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MAKE_OPAQUE(std::vector<mio::Graph<mio::SimulationNode<mio::osecirvvs::Simulation<>>, mio::MigrationEdge>>);

PYBIND11_MODULE(_simulation_osecirvvs, m)
{
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

    pymio::bind_ParameterSet<mio::osecirvvs::ParametersBase>(m, "ParametersBase");

    py::class_<mio::osecirvvs::Parameters, mio::osecirvvs::ParametersBase>(m, "Parameters")
        .def(py::init<mio::AgeGroup>())
        .def("check_constraints", &mio::osecirvvs::Parameters::check_constraints)
        .def("apply_constraints", &mio::osecirvvs::Parameters::apply_constraints);

    using SecirvvsPopulations = mio::Populations<mio::AgeGroup, mio::osecirvvs::InfectionState>;
    pymio::bind_Population(m, "Population", mio::Tag<mio::osecirvvs::Model::Populations>{});

    pymio::bind_CompartmentalModel<mio::osecirvvs::InfectionState, SecirvvsPopulations, mio::osecirvvs::Parameters>(
        m, "ModelBase");
    py::class_<mio::osecirvvs::Model, mio::CompartmentalModel<mio::osecirvvs::InfectionState, SecirvvsPopulations,
                                                              mio::osecirvvs::Parameters>>(m, "Model")
        .def(py::init<int>(), py::arg("num_agegroups"));

    pymio::bind_Simulation<mio::osecirvvs::Simulation<>>(m, "Simulation");

    m.def(
        "simulate",
        [](double t0, double tmax, double dt, const mio::osecirvvs::Model& model) {
            return mio::simulate(t0, tmax, dt, model);
        },
        "Simulates a Model from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"));

    m.def(
        "simulate_flows",
        [](double t0, double tmax, double dt, const mio::osecirvvs::Model& model) {
            return mio::simulate_flows(t0, tmax, dt, model);
        },
        "Simulates a osecirvvs model with flows from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
        py::arg("model"));

    pymio::bind_ModelNode<mio::osecirvvs::Model>(m, "ModelNode");
    pymio::bind_SimulationNode<mio::osecirvvs::Simulation<>>(m, "SimulationNode");
    pymio::bind_ModelGraph<mio::osecirvvs::Model>(m, "ModelGraph");
    pymio::bind_MigrationGraph<mio::osecirvvs::Simulation<>>(m, "MigrationGraph");
    pymio::bind_GraphSimulation<mio::Graph<mio::SimulationNode<mio::osecirvvs::Simulation<>>, mio::MigrationEdge>>(
        m, "MigrationSimulation");

    m.def(
        "draw_sample",
        [](mio::osecirvvs::Model& model) {
            return mio::osecirvvs::draw_sample(model);
        },
        py::arg("model"));

#ifdef MEMILIO_HAS_JSONCPP
    pymio::bind_write_graph<mio::osecirvvs::Model>(m);
#endif // MEMILIO_HAS_JSONCPP

    m.attr("__version__") = "dev";
}
