/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Kilian Volmer, Henrik Zunker
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
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

// Macro to skip main() in jolly_bindings.cpp
#define JOLLY_BINDINGS_SKIP_MAIN
#include "jolly_bindings.cpp"
#include "jolly_bindings2.cpp"

namespace py = pybind11;

using MobilityGraph = mio::Graph<
    mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, mio::Index<InfState>, mio::regions::Region>>,
    mio::MobilityEdgeDirected<ScalarType>>;

void bind_FarmGraph(pybind11::module_& m, std::string const& name)
{
    using G = mio::Graph<mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, mio::Index<InfState>,
                                                                        mio::regions::Region>>,
                         mio::MobilityEdgeDirected<ScalarType>>;
    pymio::bind_class<G, pymio::EnablePickling::IfAvailable>(m, name.c_str()).def(pybind11::init<>());
}

void bind_FarmGraphSimulation(pybind11::module_& m, std::string const& name)
{
    using GS = mio::FarmSimulation<
        mio::Graph<mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, mio::Index<InfState>,
                                                                  mio::regions::Region>>,
                   mio::MobilityEdgeDirected<ScalarType>>>;
    using Graph = mio::Graph<mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, mio::Index<InfState>,
                                                                            mio::regions::Region>>,
                             mio::MobilityEdgeDirected<ScalarType>>;
    pymio::bind_class<GS, pymio::EnablePickling::Never>(m, name.c_str())
        .def(pybind11::init([](Graph& graph, double t0, double dt, u_int seed) {
                 auto rng = mio::RandomNumberGenerator();
                 rng.seed({seed});
                 return std::make_unique<GS>(mio::make_farm_sim(t0, dt, std::move(graph), rng));
             }),
             pybind11::arg("graph"), pybind11::arg("t0") = 0.0, pybind11::arg("dt") = 1.0, pybind11::arg("seed") = 0);
}

PYBIND11_MODULE(_simulation_jolly, m)
{
    m.def("simulate", &simulate, "Simulates the jolly model", py::arg("farm_file"), py::arg("tmax") = 40,
          py::arg("dt") = 1.0, py::arg("suspicion_threshold") = 0.05, py::arg("sensitivity") = 0.95,
          py::arg("h0") = 0.001, py::arg("r0") = 1000, py::arg("alpha") = 2, py::arg("infection_baseline") = 0.0001,
          py::arg("culling_factor") = 0.99, py::arg("A0_SEI") = 0.2, py::arg("A0_EI") = 0.2, py::arg("A0_ID") = 0.2,
          py::arg("A0_DeathRate") = 0.000128, py::arg("A1_SEI") = 0.2, py::arg("A1_EI") = 0.2, py::arg("A1_ID") = 0.2,
          py::arg("A1_DeathRate") = 0.000128, py::arg("A2_SEI") = 0.2, py::arg("A2_EI") = 0.2, py::arg("A2_ID") = 0.2,
          py::arg("A2_DeathRate") = 0.000588, py::arg("A3_SEI") = 0.2, py::arg("A3_EI") = 0.2, py::arg("A3_ID") = 0.2,
          py::arg("A3_DeathRate") = 0.000588, py::arg("A4_SEI") = 0.2, py::arg("A4_EI") = 0.2, py::arg("A4_ID") = 0.2,
          py::arg("A4_DeathRate") = 0.000588, py::arg("foi_inner_factor0") = 1.0, py::arg("foi_outer_factor0") = 1.0,
          py::arg("foi_inner_factor1") = 1.0, py::arg("foi_outer_factor1") = 1.0, py::arg("foi_inner_factor2") = 1.0,
          py::arg("foi_outer_factor2") = 1.0, py::arg("foi_inner_factor3") = 1.0, py::arg("foi_outer_factor3") = 1.0,
          py::arg("foi_inner_factor4") = 1.0, py::arg("foi_outer_factor4") = 1.0, py::arg("damping0") = 1.0,
          py::arg("damping1") = 1.0, py::arg("damping2") = 1.0, py::arg("damping3") = 1.0, py::arg("damping4") = 1.0,
          py::arg("first_infection_day") = 0, py::arg("second_infection_day") = 2, py::arg("third_infection_day") = 2,
          py::arg("seed") = 42);

    m.def("simulate_with_init", &simulate_with_init, "Simulates the jolly model with initialization from farm file",
          py::arg("farm_file"), py::arg("tmax"), py::arg("dt"), py::arg("suspicion_threshold"), py::arg("sensitivity"),
          py::arg("h0"), py::arg("r0"), py::arg("alpha"), py::arg("infection_baseline"), py::arg("culling_factor"),
          py::arg("A0_SEI"), py::arg("A0_EI"), py::arg("A0_ID"), py::arg("A0_DeathRate"), py::arg("A1_SEI"),
          py::arg("A1_EI"), py::arg("A1_ID"), py::arg("A1_DeathRate"), py::arg("A2_SEI"), py::arg("A2_EI"),
          py::arg("A2_ID"), py::arg("A2_DeathRate"), py::arg("A3_SEI"), py::arg("A3_EI"), py::arg("A3_ID"),
          py::arg("A3_DeathRate"), py::arg("A4_SEI"), py::arg("A4_EI"), py::arg("A4_ID"), py::arg("A4_DeathRate"),
          py::arg("foi_inner_factor0") = 1.0, py::arg("foi_outer_factor0"), py::arg("foi_inner_factor1"),
          py::arg("foi_outer_factor1"), py::arg("foi_inner_factor2"), py::arg("foi_outer_factor2"),
          py::arg("foi_inner_factor3"), py::arg("foi_outer_factor3"), py::arg("foi_inner_factor4"),
          py::arg("foi_outer_factor4"), py::arg("damping0"), py::arg("damping1"), py::arg("damping2"),
          py::arg("damping3"), py::arg("damping4"), py::arg("first_infection_day"), py::arg("second_infection_day"),
          py::arg("third_infection_day"), py::arg("seed"), pybind11::return_value_policy::reference_internal);

    m.def("simulate_continued", &simulate_continued, "Simulates the jolly model continued from other simulation",
          py::arg("sim_to_copy"), py::arg("t_start"), py::arg("tmax"), py::arg("dt"), py::arg("suspicion_threshold"),
          py::arg("sensitivity"), py::arg("h0"), py::arg("r0"), py::arg("alpha"), py::arg("infection_baseline"),
          py::arg("culling_factor"), py::arg("A0_SEI"), py::arg("A0_EI"), py::arg("A0_ID"), py::arg("A0_DeathRate"),
          py::arg("A1_SEI"), py::arg("A1_EI"), py::arg("A1_ID"), py::arg("A1_DeathRate"), py::arg("A2_SEI"),
          py::arg("A2_EI"), py::arg("A2_ID"), py::arg("A2_DeathRate"), py::arg("A3_SEI"), py::arg("A3_EI"),
          py::arg("A3_ID"), py::arg("A3_DeathRate"), py::arg("A4_SEI"), py::arg("A4_EI"), py::arg("A4_ID"),
          py::arg("A4_DeathRate"), py::arg("foi_inner_factor0") = 1.0, py::arg("foi_outer_factor0"),
          py::arg("foi_inner_factor1"), py::arg("foi_outer_factor1"), py::arg("foi_inner_factor2"),
          py::arg("foi_outer_factor2"), py::arg("foi_inner_factor3"), py::arg("foi_outer_factor3"),
          py::arg("foi_inner_factor4"), py::arg("foi_outer_factor4"), py::arg("damping0"), py::arg("damping1"),
          py::arg("damping2"), py::arg("damping3"), py::arg("damping4"), py::arg("first_infection_day"),
          py::arg("second_infection_day"), py::arg("third_infection_day"), py::arg("seed"),
          pybind11::return_value_policy::reference_internal);

    m.def(
        "get_result",
        [](mio::FarmSimulation<
               mio::Graph<mio::FarmNode<ScalarType, mio::smm::Simulation<ScalarType, InfState, mio::Index<InfState>,
                                                                         mio::regions::Region>>,
                          mio::MobilityEdgeDirected<ScalarType>>>& graph_sim,
           double tmax) {
            return graph_sim.get_confirmation_dates(tmax);
        },
        py::return_value_policy::reference_internal);

    bind_FarmGraph(m, "FarmGraph");
    bind_FarmGraphSimulation(m, "FarmSimulation");

    m.attr("__version__") = "dev";
}
