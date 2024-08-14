/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "simulation.h"
#include "models/osir.h"
#include "models/oseir.h"
#include "models/osecir.h"
#include "models/osecirvvs.h"
#include "models/abm.h"

#include "pybind11/pybind11.h"

PYBIND11_MODULE(simulation, m_simulation)
{
   pymio::bind_simulation(m_simulation);

   auto m_osir = m_simulation.def_submodule("osir");
   pymio::bind_osir(m_osir);

   auto m_oseir = m_simulation.def_submodule("oseir");
   pymio::bind_oseir(m_oseir);

   auto m_osecir = m_simulation.def_submodule("osecir");
   pymio::bind_osecir(m_osecir);

   auto m_osecirvvs = m_simulation.def_submodule("osecirvvs");
   pymio::bind_osecirvvs(m_osecirvvs);

   auto m_abm = m_simulation.def_submodule("abm");
   pymio::bind_abm(m_abm);
}