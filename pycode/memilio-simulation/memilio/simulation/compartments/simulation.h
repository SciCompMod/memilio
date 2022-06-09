#ifndef PYMIO_SIMULATION_H
#define PYMIO_SIMULATION_H

#include "memilio/compartments/simulation.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{

/*
 * @brief bind Simulation for any number model
 */
template <class Simulation>
void bind_Simulation(pybind11::module& m, std::string const& name)
{
    pybind11::class_<Simulation>(m, name.c_str())
        .def(pybind11::init<const typename Simulation::Model&, double, double>(), pybind11::arg("model"), pybind11::arg("t0") = 0,
             pybind11::arg("dt") = 0.1)
        .def_property_readonly("result", pybind11::overload_cast<>(&Simulation::get_result, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def_property_readonly("model", pybind11::overload_cast<>(&Simulation::get_model, pybind11::const_),
                               pybind11::return_value_policy::reference_internal)
        .def("advance", &Simulation::advance, pybind11::arg("tmax"));
}

} // namespace pymio

#endif //PYMIO_SIMULATION_H