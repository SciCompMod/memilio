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
#include "utils/index.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "compartments/simulation.h"
#include "compartments/compartmentalmodel.h"
#include "epidemiology/populations.h"
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"

namespace py = pybind11;

namespace pymio{
//specialization of pretty_name
template <>
std::string pretty_name<mio::oseir::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MODULE(_simulation_oseir, m)
{
pymio::iterable_enum<mio::oseir::InfectionState>(m, "InfectionState")
    .value("Susceptible", mio::oseir::InfectionState::Susceptible)
    .value("Exposed", mio::oseir::InfectionState::Exposed)
    .value("Infected", mio::oseir::InfectionState::Infected)
    .value("Recovered", mio::oseir::InfectionState::Recovered);

pymio::bind_Index<mio::oseir::InfectionState>(m, "Index_InfectionState");
pymio::bind_CustomIndexArray<mio::UncertainValue, mio::oseir::InfectionState>(m, "PopulationArray");

pymio::bind_ParameterSet<mio::oseir::Parameters>(m, "Parameters");

pymio::bind_Population<mio::oseir::InfectionState>(m, "Population");

using Populations = mio::Populations<mio::oseir::InfectionState>;
pymio::bind_CompartmentalModel<mio::oseir::InfectionState, Populations, mio::oseir::Parameters>(m, "ModelBase");
py::class_<mio::oseir::Model, mio::CompartmentalModel<mio::oseir::InfectionState, Populations, mio::oseir::Parameters>>(m, "Model")
   .def(py::init<>());

m.def(
"simulate",
    [](double t0, double tmax, double dt, const mio::oseir::Model& model) {
        return mio::simulate(t0, tmax, dt, model);
    },
    "Simulates a oseir from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"));

m.attr("__version__") = "dev";
}
