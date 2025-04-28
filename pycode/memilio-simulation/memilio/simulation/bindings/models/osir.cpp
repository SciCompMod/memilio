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
#include "ode_sir/model.h"
#include "ode_sir/infection_state.h"
#include "memilio/data/analyze_result.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{
//specialization of pretty_name
template <>
inline std::string pretty_name<mio::osir::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MODULE(_simulation_osir, m)
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

    pymio::iterable_enum<mio::osir::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::osir::InfectionState::Susceptible)
        .value("Infected", mio::osir::InfectionState::Infected)
        .value("Recovered", mio::osir::InfectionState::Recovered);

    pymio::bind_ParameterSet<mio::osir::ParametersBase<double>, pymio::EnablePickling::Required>(m, "ParametersBase");

    pymio::bind_class<mio::osir::Parameters<double>, pymio::EnablePickling::Required,
                      mio::osir::ParametersBase<double>>(m, "Parameters")
        .def(py::init<mio::AgeGroup>())
        .def("check_constraints", &mio::osir::Parameters<double>::check_constraints);

    using Populations = mio::Populations<double, mio::AgeGroup, mio::osir::InfectionState>;
    pymio::bind_Population(m, "Populations", mio::Tag<mio::osir::Model<double>::Populations>{});
    pymio::bind_CompartmentalModel<mio::osir::InfectionState, Populations, mio::osir::Parameters<double>,
                                   pymio::EnablePickling::Never>(m, "ModelBase");
    pymio::bind_class<
        mio::osir::Model<double>, pymio::EnablePickling::Required,
        mio::CompartmentalModel<double, mio::osir::InfectionState, Populations, mio::osir::Parameters<double>>>(
        m, "Model")
        .def(py::init<int>(), py::arg("num_agegroups"));

    pymio::bind_Simulation<mio::Simulation<double, mio::osir::Model<double>>>(m, "Simulation");
    
    m.def(
        "simulate", &mio::simulate<double, mio::osir::Model<double>>, "Simulates an ODE SIR model from t0 to tmax.", 
        py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none());

    m.attr("__version__") = "dev";
}
