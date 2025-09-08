/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Maximilian Betz
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
#include "compartments/compartmental_model.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "epidemiology/populations.h"

//Includes from MEmilio
#include "sde_sir/model.h"
#include "sde_sir/infection_state.h"
#include "memilio/compartments/stochastic_simulation.h"
#include "memilio/compartments/stochastic_model.h"
#include "memilio/data/analyze_result.h"

#include "pybind11/pybind11.h"
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;

namespace pymio
{
//specialization of pretty_name
template <>
inline std::string pretty_name<mio::ssir::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MODULE(_simulation_ssir, m)
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

    pymio::iterable_enum<mio::ssir::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::ssir::InfectionState::Susceptible)
        .value("Infected", mio::ssir::InfectionState::Infected)
        .value("Recovered", mio::ssir::InfectionState::Recovered);

    pymio::bind_ParameterSet<mio::ssir::ParametersBase, pymio::EnablePickling::Never>(m, "ParametersBase");

    pymio::bind_class<mio::ssir::Parameters, pymio::EnablePickling::Required, mio::ssir::ParametersBase>(m,
                                                                                                         "Parameters")
        .def(py::init<>())
        .def("check_constraints", &mio::ssir::Parameters::check_constraints)
        .def("apply_constraints", &mio::ssir::Parameters::apply_constraints);

    using Populations = mio::Populations<double, mio::ssir::InfectionState>;
    pymio::bind_Population(m, "Populations", mio::Tag<mio::ssir::Model::Populations>{});
    pymio::bind_CompartmentalModel<mio::ssir::InfectionState, Populations, mio::ssir::Parameters,
                                   pymio::EnablePickling::Never>(m, "ModelBase");
    // pymio::bind_StochasticModel<double, mio::ssir::InfectionState, Populations, mio::ssir::Parameters,
    //                                pymio::EnablePickling::Never>(m, "ModelBase");
    pymio::bind_class<mio::ssir::Model, pymio::EnablePickling::Never,
                      mio::CompartmentalModel<double, mio::ssir::InfectionState, Populations, mio::ssir::Parameters>>(
        m, "Model")
        .def(py::init<>());

    // pymio::bind_Simulation<mio::Simulation<double, mio::ssir::Model>>(m, "Simulation");

    m.def(
        "simulate_stochastic",
        [](double t0, double tmax, double dt, mio::ssir::Model const& model) {
            return mio::simulate_stochastic<double, mio::ssir::Model>(t0, tmax, dt, model);
        },
        "Simulates an SDE SIR model from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"));

    m.attr("__version__") = "dev";
}
