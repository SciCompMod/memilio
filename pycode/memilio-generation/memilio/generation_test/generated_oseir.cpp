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

#include "templates.h"
#include "ode_seir/model.h"
#include "Eigen/Core"
#include "pybind11/stl_bind.h"
#include <vector>


namespace py = pybind11;

namespace pymio{
//specialization of pretty_name
template <>
std::string pretty_name<mio::oseir::InfectionState>()
{
   return "InfectionState";
}

} // namespace pymio


PYBIND11_MAKE_OPAQUE(std::vector<mio::Graph<mio::SimulationNode<${simulation_class}>, mio::MigrationEdge>>);

PYBIND11_MODULE(_simulation_generated_oseir, m)
{
    pymio::iterable_enum<mio::oseir::InfectionState>(m, "InfectionState")
	    .value("InfectionState::Susceptible", mio::oseir::InfectionState::Susceptible)
	    .value("InfectionState::Exposed", mio::oseir::InfectionState::Exposed)
	    .value("InfectionState::Infected", mio::oseir::InfectionState::Infected)
	    .value("InfectionState::Recovered", mio::oseir::InfectionState::Recovered)
	    .value("InfectionState::Count", mio::oseir::InfectionState::Count);



    None
    pymio::bind_ParameterSet<mio::oseir::Parameters>(m, "Parameters");

    None

    using Populations = mio::Populations<mio::oseir::InfectionState>;
    pymio::bind_Population(m, "Population", mio::Tag<mio::oseir::Model::Populations>{});
    
    pymio::bind_CompartmentalModel<mio::Populations<mio::oseir::InfectionState>, mio::ParameterSet<mio::oseir::InfectionProbabilityFromContact, mio::oseir::LatentTime, mio::oseir::InfectiousTime, mio::oseir::ContactPatterns>>(m, "ModelBase");
    py::class_<mio::oseir::Model, mio::CompartmentalModel<mio::Populations<mio::oseir::InfectionState>, mio::ParameterSet<mio::oseir::InfectionProbabilityFromContact, mio::oseir::LatentTime, mio::oseir::InfectiousTime, mio::oseir::ContactPatterns>>>(m, "Model")
        .def(py::init<>());


    

    m.def(
        "simulate",
        [](double t0, double tmax, double dt, const mio::oseir::Model& model) {
            return mio::simulate(t0, tmax, dt, model);
        },
        "Simulates a Model from t0 to tmax.", py::arg("t0"), py::arg("tmax"), py::arg("dt"), py::arg("model"));
    
    

    m.attr("__version__") = "dev";
}
