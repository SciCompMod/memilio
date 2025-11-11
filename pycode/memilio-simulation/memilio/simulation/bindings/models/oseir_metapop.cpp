/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "epidemiology/age_group.h"
#include "epidemiology/populations.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "utils/index.h"

// Includes from MEmilio
#include "models/ode_seir_metapop/model.h"
#include "models/ode_seir_metapop/infection_state.h"
#include "models/ode_seir_metapop/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"

#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"

namespace py = pybind11;

namespace pymio
{
// specialization of pretty_name
template <>
inline std::string pretty_name<mio::oseirmetapop::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

PYBIND11_MODULE(_simulation_oseir_metapop, m)
{
    pymio::iterable_enum<mio::oseirmetapop::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::oseirmetapop::InfectionState::Susceptible)
        .value("Exposed", mio::oseirmetapop::InfectionState::Exposed)
        .value("Infected", mio::oseirmetapop::InfectionState::Infected)
        .value("Recovered", mio::oseirmetapop::InfectionState::Recovered);

    pymio::bind_ParameterSet<mio::oseirmetapop::ParametersBase<double>, pymio::EnablePickling::Required>(
        m, "ParametersBase");

    pymio::bind_class<mio::oseirmetapop::Parameters<double>, pymio::EnablePickling::Never,
                      mio::oseirmetapop::ParametersBase<double>>(m, "Parameters")
        .def(py::init<mio::regions::Region, mio::AgeGroup>(), py::arg("num_regions"), py::arg("num_agegroups"))
        .def("check_constraints", &mio::oseirmetapop::Parameters<double>::check_constraints)
        .def("apply_constraints", &mio::oseirmetapop::Parameters<double>::apply_constraints)
        .def_property_readonly("num_regions", &mio::oseirmetapop::Parameters<double>::get_num_regions)
        .def_property_readonly("num_agegroups", &mio::oseirmetapop::Parameters<double>::get_num_agegroups);

    using Populations =
        mio::Populations<double, mio::regions::Region, mio::AgeGroup, mio::oseirmetapop::InfectionState>;

    pymio::bind_Population(m, "Populations", mio::Tag<mio::oseirmetapop::Model<double>::Populations>{});
    pymio::bind_CompartmentalModel<mio::oseirmetapop::InfectionState, Populations,
                                   mio::oseirmetapop::Parameters<double>, pymio::EnablePickling::Never>(m, "ModelBase");

    using Model = mio::oseirmetapop::Model<double>;
    pymio::bind_class<Model, pymio::EnablePickling::Never,
                      mio::CompartmentalModel<double, mio::oseirmetapop::InfectionState, Populations,
                                              mio::oseirmetapop::Parameters<double>>>(m, "Model")
        .def(py::init<int, int>(), py::arg("num_regions"), py::arg("num_agegroups"))
        .def("set_commuting_strengths", py::overload_cast<const Eigen::MatrixXd&>(&Model::set_commuting_strengths),
             py::arg("commuting_strengths"))
        .def("set_commuting_strengths_identity", py::overload_cast<>(&Model::set_commuting_strengths));

    pymio::bind_Simulation<mio::Simulation<double, Model>>(m, "Simulation");

    m.def("simulate", &mio::simulate<double, Model>, "Simulates an ODE SEIR metapopulation model from t0 to tmax.",
          py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"), py::arg("integrator") = py::none());

    m.attr("__version__") = "dev";
}
