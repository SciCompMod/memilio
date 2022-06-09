#ifndef PYMIO_GRAPH_SIMULATION_H
#define PYMIO_GRAPH_SIMULATION_H

#include "memilio/mobility/mobility.h"
#include "memilio/mobility/graph_simulation.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{

/*
 * @brief bind GraphSimulation for any node and edge type
 */
template <class Graph>
void bind_GraphSimulation(pybind11::module& m, std::string const& name)
{
    using GS = mio::GraphSimulation<Graph>;
    pybind11::class_<GS>(m, name.c_str())
        .def(pybind11::init([](Graph& graph, double t0, double dt) {
                 return std::make_unique<GS>(mio::make_migration_sim(t0, dt, std::move(graph)));
             }),
             pybind11::arg("graph"), pybind11::arg("t0") = 0.0, pybind11::arg("dt") = 1.0)
        .def_property_readonly(
            "graph", [](GS& self) -> Graph& { return self.get_graph(); },
            pybind11::return_value_policy::reference_internal)
        .def_property_readonly("t", &GS::get_t)
        .def("advance", &GS::advance, pybind11::arg("tmax"));
}


} // namespace pymio

#endif //PYMIO_GRAPH_SIMULATION_H