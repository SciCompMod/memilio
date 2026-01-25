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
#include "sde_seirvv/model.h"
#include "sde_seirvv/infection_state.h"
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
inline std::string pretty_name<mio::sseirvv::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MODULE(_simulation_sseirvv, m)
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

    pymio::iterable_enum<mio::sseirvv::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::sseirvv::InfectionState::Susceptible)
        .value("ExposedV1", mio::sseirvv::InfectionState::ExposedV1)
        .value("InfectedV1", mio::sseirvv::InfectionState::InfectedV1)
        .value("RecoveredV1", mio::sseirvv::InfectionState::RecoveredV1)
        .value("ExposedV2", mio::sseirvv::InfectionState::ExposedV2)
        .value("InfectedV2", mio::sseirvv::InfectionState::InfectedV2)
        .value("RecoveredV2", mio::sseirvv::InfectionState::RecoveredV2)
        .value("ExposedV1V2", mio::sseirvv::InfectionState::ExposedV1V2)
        .value("InfectedV1V2", mio::sseirvv::InfectionState::InfectedV1V2)
        .value("RecoveredV1V2", mio::sseirvv::InfectionState::RecoveredV1V2);

    pymio::bind_ParameterSet<mio::sseirvv::ParametersBase<double>, pymio::EnablePickling::Never>(m, "ParametersBase");

    pymio::bind_class<mio::sseirvv::Parameters<double>, pymio::EnablePickling::Required,
                      mio::sseirvv::ParametersBase<double>>(m, "Parameters")
        .def(py::init<>())
        .def("check_constraints", &mio::sseirvv::Parameters<double>::check_constraints)
        .def("apply_constraints", &mio::sseirvv::Parameters<double>::apply_constraints);

    using Populations = mio::Populations<double, mio::sseirvv::InfectionState>;
    pymio::bind_Population(m, "Populations", mio::Tag<mio::sseirvv::Model<double>::Populations>{});
    pymio::bind_CompartmentalModel<mio::sseirvv::InfectionState, Populations, mio::sseirvv::Parameters<double>,
                                   pymio::EnablePickling::Never>(m, "ModelBase");
    pymio::bind_class<
        mio::sseirvv::Model<double>, pymio::EnablePickling::Never,
        mio::CompartmentalModel<double, mio::sseirvv::InfectionState, Populations, mio::sseirvv::Parameters<double>>>(
        m, "Model")
        .def(py::init<>());

    m.def(
        "simulate_stochastic",
        [](double t0, double tmax, double dt, mio::sseirvv::Model<double> const& model) {
            return mio::simulate_stochastic<double, mio::sseirvv::Model<double>>(t0, tmax, dt, model);
        },
        "Simulates an SDE SEIRVV model from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
        py::arg("model"));

    m.attr("__version__") = "dev";
}
