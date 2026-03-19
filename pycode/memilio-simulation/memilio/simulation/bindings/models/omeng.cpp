/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
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

// Includes from pymio
#include "pybind_util.h"
#include "compartments/simulation.h"
#include "compartments/flow_simulation.h"
#include "compartments/compartmental_model.h"
#include "epidemiology/age_group.h"
#include "epidemiology/populations.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "utils/index.h"
#include "data/analyze_result.h"

// Includes from MEmilio
#include "ode_meningitis/model.h"
#include "ode_meningitis/infection_state.h"
#include "ode_meningitis/parameters.h"
#include "memilio/data/analyze_result.h"

#include "pybind11/pybind11.h"
#include "Eigen/Core"

namespace py = pybind11;

namespace pymio
{
// specialization of pretty_name
template <>
inline std::string pretty_name<mio::omeng::InfectionState>()
{
    return "InfectionState";
}
} // namespace pymio

PYBIND11_MODULE(_simulation_omeng, m)
{
    // interpolation helpers (interpolate_simulation_result)
    pymio::bind_interpolate_result_methods(m);

    // InfectionState enum
    pymio::iterable_enum<mio::omeng::InfectionState>(m, "InfectionState")
        .value("Incoming", mio::omeng::InfectionState::Incoming)
        .value("SusceptibleLow", mio::omeng::InfectionState::SusceptibleLow)
        .value("SusceptibleHigh", mio::omeng::InfectionState::SusceptibleHigh)
        .value("Carrier", mio::omeng::InfectionState::Carrier)
        .value("Infected", mio::omeng::InfectionState::Infected)
        .value("Recovered", mio::omeng::InfectionState::Recovered)
        .value("Dead", mio::omeng::InfectionState::Dead)
        .value("DeadNatural", mio::omeng::InfectionState::DeadNatural)
        .value("Count", mio::omeng::InfectionState::Count);

    // Parameters
    pymio::bind_ParameterSet<mio::omeng::ParametersBase<double>, pymio::EnablePickling::Required>(m, "ParametersBase");

    pymio::bind_class<mio::omeng::Parameters<double>, pymio::EnablePickling::Required,
                      mio::omeng::ParametersBase<double>>(m, "Parameters")
        .def(py::init<mio::AgeGroup>())
        .def("check_constraints", &mio::omeng::Parameters<double>::check_constraints)
        .def("apply_constraints", &mio::omeng::Parameters<double>::apply_constraints);

    // Populations and Model
    using Populations = mio::Populations<double, mio::AgeGroup, mio::omeng::InfectionState>;
    pymio::bind_Population(m, "Populations", mio::Tag<mio::omeng::Model<double>::Populations>{});

    // Bind the CompartmentalModel base (FlowModel inherits from CompartmentalModel)
    pymio::bind_CompartmentalModel<mio::omeng::InfectionState, Populations, mio::omeng::Parameters<double>,
                                   pymio::EnablePickling::Never>(m, "ModelBase");

    pymio::bind_class<
        mio::omeng::Model<double>, pymio::EnablePickling::Required,
        mio::CompartmentalModel<double, mio::omeng::InfectionState, Populations, mio::omeng::Parameters<double>>>(
        m, "Model")
        .def(py::init<int>(), py::arg("num_agegroups"));

    // Simulation (with NPI support via overridden advance())
    pymio::bind_Simulation<mio::omeng::Simulation<double>>(m, "Simulation");

    // FlowSimulation (also records all flows at each time step)
    pymio::bind_Flow_Simulation<mio::FlowSimulation<double, mio::omeng::Model<double>>>(m, "FlowSimulation");

    // simulate: uses the model-specific Simulation class (includes dynamic NPI checks in advance())
    m.def("simulate", &mio::omeng::simulate<double>, "Simulates the ODE meningitis model from t0 to tmax.",
          py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none());

    // simulate_flows: additionally records all inter-compartment flows at each time step.
    // Returns a pair (compartments TimeSeries, flows TimeSeries).
    m.def("simulate_flows", &mio::omeng::simulate_flows<double>,
          "Simulates the ODE meningitis model with flows from t0 to tmax.", py::arg("t0"), py::arg("tmax"),
          py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none());

    m.attr("__version__") = "dev";
}
