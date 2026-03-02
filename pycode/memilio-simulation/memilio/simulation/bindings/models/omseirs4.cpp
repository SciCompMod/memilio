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
#include "compartments/compartmental_model.h"
#include "epidemiology/populations.h"
#include "utils/parameter_set.h"
#include "utils/index.h"
#include "data/analyze_result.h"

// Includes from MEmilio
#include "ode_mseirs4/model.h"
#include "ode_mseirs4/infection_state.h"
#include "ode_mseirs4/parameters.h"
#include "memilio/data/analyze_result.h"

#include "pybind11/pybind11.h"
#include "Eigen/Core"
#include <vector>

namespace py = pybind11;

namespace pymio
{
// specialization of pretty_name
template <>
inline std::string pretty_name<mio::omseirs4::InfectionState>()
{
    return "InfectionState";
}
} // namespace pymio

PYBIND11_MODULE(_simulation_omseirs4, m)
{
    // interpolation helpers
    pymio::bind_interpolate_result_methods(m);

    // InfectionState enum
    pymio::iterable_enum<mio::omseirs4::InfectionState>(m, "InfectionState")
        .value("MaternalImmune", mio::omseirs4::InfectionState::MaternalImmune)
        .value("S1", mio::omseirs4::InfectionState::S1)
        .value("S2", mio::omseirs4::InfectionState::S2)
        .value("S3", mio::omseirs4::InfectionState::S3)
        .value("S4", mio::omseirs4::InfectionState::S4)
        .value("E1", mio::omseirs4::InfectionState::E1)
        .value("E2", mio::omseirs4::InfectionState::E2)
        .value("E3", mio::omseirs4::InfectionState::E3)
        .value("E4", mio::omseirs4::InfectionState::E4)
        .value("I1", mio::omseirs4::InfectionState::I1)
        .value("I2", mio::omseirs4::InfectionState::I2)
        .value("I3", mio::omseirs4::InfectionState::I3)
        .value("I4", mio::omseirs4::InfectionState::I4)
        .value("R1", mio::omseirs4::InfectionState::R1)
        .value("R2", mio::omseirs4::InfectionState::R2)
        .value("R3", mio::omseirs4::InfectionState::R3)
        .value("R4", mio::omseirs4::InfectionState::R4);

    // Parameters
    pymio::bind_ParameterSet<mio::omseirs4::ParametersBase<double>, pymio::EnablePickling::Required>(m,
                                                                                                     "ParametersBase");

    pymio::bind_class<mio::omseirs4::Parameters<double>, pymio::EnablePickling::Required,
                      mio::omseirs4::ParametersBase<double>>(m, "Parameters")
        .def(py::init<>())
        .def("check_constraints", &mio::omseirs4::Parameters<double>::check_constraints)
        .def("apply_constraints", &mio::omseirs4::Parameters<double>::apply_constraints);

    // Populations and Model
    using Populations = mio::Populations<double, mio::omseirs4::InfectionState>;
    pymio::bind_Population(m, "Populations", mio::Tag<Populations>{});
    pymio::bind_CompartmentalModel<mio::omseirs4::InfectionState, Populations, mio::omseirs4::Parameters<double>,
                                   pymio::EnablePickling::Never>(m, "ModelBase");
    pymio::bind_class<
        mio::omseirs4::Model<double>, pymio::EnablePickling::Required,
        mio::CompartmentalModel<double, mio::omseirs4::InfectionState, Populations, mio::omseirs4::Parameters<double>>>(
        m, "Model")
        .def(py::init<>());

    // Simulation and simulate()
    using Sim = mio::Simulation<double, mio::omseirs4::Model<double>>;
    pymio::bind_Simulation<Sim>(m, "Simulation");

    m.def("simulate", &mio::simulate<double, mio::omseirs4::Model<double>>,
          "Simulates an ODE MSEIRS4 model from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"),
          py::arg("model"), py::arg("integrator") = py::none());

    m.attr("__version__") = "dev";
}
