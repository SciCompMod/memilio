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
#include "utils/index.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "compartments/simulation.h"
#include "compartments/flow_simulation.h"
#include "compartments/compartmentalmodel.h"
#include "epidemiology/age_group.h"
#include "epidemiology/populations.h"

//Includes from MEmilio
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "memilio/data/analyze_result.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{
//specialization of pretty_name
template <>
inline std::string pretty_name<mio::oseir::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MODULE(_simulation_oseir, m)
{
    m.def("interpolate_simulation_result",
          static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&, const double)>(
              &mio::interpolate_simulation_result),
          py::arg("ts"), py::arg("abs_tol") = 1e-14);

    m.def("interpolate_simulation_result",
          static_cast<mio::TimeSeries<double> (*)(const mio::TimeSeries<double>&, const std::vector<double>&)>(
              &mio::interpolate_simulation_result),
          py::arg("ts"), py::arg("interpolation_times"));

    m.def("interpolate_ensemble_results", &mio::interpolate_ensemble_results<mio::TimeSeries<double>>);

    pymio::iterable_enum<mio::oseir::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::oseir::InfectionState::Susceptible)
        .value("Exposed", mio::oseir::InfectionState::Exposed)
        .value("Infected", mio::oseir::InfectionState::Infected)
        .value("Recovered", mio::oseir::InfectionState::Recovered);

    pymio::bind_ParameterSet<mio::oseir::ParametersBase<double>, pymio::EnablePickling::Required>(m, "ParametersBase");

    pymio::bind_class<mio::oseir::Parameters<double>, pymio::EnablePickling::Required,
                      mio::oseir::ParametersBase<double>>(m, "Parameters", py::module_local{})

        .def(py::init<mio::AgeGroup>())
        .def("check_constraints", &mio::oseir::Parameters<double>::check_constraints);

    using Populations = mio::Populations<double, mio::AgeGroup, mio::oseir::InfectionState>;
    pymio::bind_Population(m, "Populations", mio::Tag<mio::oseir::Model<double>::Populations>{});
    pymio::bind_CompartmentalModel<mio::oseir::InfectionState, Populations, mio::oseir::Parameters<double>,
                                   pymio::EnablePickling::Never>(m, "ModelBase");
    pymio::bind_class<
        mio::oseir::Model<double>, pymio::EnablePickling::Required,
        mio::CompartmentalModel<double, mio::oseir::InfectionState, Populations, mio::oseir::Parameters<double>>>(
        m, "Model")
        .def(py::init<int>(), py::arg("num_agegroups"));

    pymio::bind_Simulation<mio::Simulation<double, mio::oseir::Model<double>>>(m, "Simulation");
    pymio::bind_Flow_Simulation<mio::FlowSimulation<double, mio::oseir::Model<double>>>(m, "FlowSimulation");

    m.def(
        "simulate", &mio::simulate<double, mio::oseir::Model<double>>, "Simulates an ODE SEIR from t0 to tmax.", 
        py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none());

    m.def(
        "simulate_flows", &mio::simulate_flows<double, mio::oseir::Model<double>>, "Simulates an ODE SEIR with flows from t0 to tmax.", 
        py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none());

    m.attr("__version__") = "dev";
}
