/*
* Copyright (C) 2020-2026 MEmilio
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
#include "sde_sirs/model.h"
#include "sde_sirs/infection_state.h"
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
inline std::string pretty_name<mio::ssirs::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MODULE(_simulation_ssirs, m)
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

    pymio::iterable_enum<mio::ssirs::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::ssirs::InfectionState::Susceptible)
        .value("Infected", mio::ssirs::InfectionState::Infected)
        .value("Recovered", mio::ssirs::InfectionState::Recovered);

    pymio::bind_ParameterSet<mio::ssirs::ParametersBase<double>, pymio::EnablePickling::Never>(m, "ParametersBase");

    pymio::bind_class<mio::ssirs::Parameters<double>, pymio::EnablePickling::Required,
                      mio::ssirs::ParametersBase<double>>(m, "Parameters")
        .def(py::init<>())
        .def("check_constraints", &mio::ssirs::Parameters<double>::check_constraints)
        .def("apply_constraints", &mio::ssirs::Parameters<double>::apply_constraints);

    using Populations = mio::Populations<double, mio::ssirs::InfectionState>;
    pymio::bind_Population(m, "Populations", mio::Tag<mio::ssirs::Model<double>::Populations>{});
    pymio::bind_CompartmentalModel<mio::ssirs::InfectionState, Populations, mio::ssirs::Parameters<double>,
                                   pymio::EnablePickling::Never>(m, "ModelBase");
    pymio::bind_class<
        mio::ssirs::Model<double>, pymio::EnablePickling::Never,
        mio::CompartmentalModel<double, mio::ssirs::InfectionState, Populations, mio::ssirs::Parameters<double>>>(
        m, "Model")
        .def(py::init<>());

    m.def(
        "simulate_stochastic",
        [](double t0, double tmax, double dt, mio::ssirs::Model<double> const& model) {
            return mio::simulate_stochastic<double, mio::ssirs::Model<double>>(t0, tmax, dt, model);
        },
        "Simulates an SDE SIRS model from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
        py::arg("model"));

    m.attr("__version__") = "dev";
}
