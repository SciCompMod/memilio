#include "pybind11/pybind11.h"

namespace py = pybind11;

// special submodule
void bind_simulation(py::module &);
void bind_osir(py::module &);
void bind_oseir(py::module &);
void bind_osecir(py::module &);
void bind_osecirvvs(py::module &);
void bind_abm(py::module &);

PYBIND11_MODULE(simulation, m_simulation)
{
   bind_simulation(m_simulation);

   auto m_osir = m_simulation.def_submodule("osir");
   bind_osir(m_osir);

   auto m_oseir = m_simulation.def_submodule("oseir");
   bind_oseir(m_oseir);

   auto m_osecir = m_simulation.def_submodule("osecir");
   bind_osecir(m_osecir);

   auto m_osecirvvs = m_simulation.def_submodule("osecirvvs");
   bind_osecirvvs(m_osecirvvs);

   auto m_abm = m_simulation.def_submodule("abm");
   bind_abm(m_abm);
}