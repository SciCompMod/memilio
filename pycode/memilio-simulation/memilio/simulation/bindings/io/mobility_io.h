/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef PYMIO_MOBILITY_IO_H
#define PYMIO_MOBILITY_IO_H

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/io/mobility_io.h"

#include "pybind11/pybind11.h"
#include <cstddef>

namespace pymio
{

template <class Model>
void bind_write_graph(pybind11::module_& m)
{
    m.def("write_graph",
          [&](const mio::Graph<Model, mio::MobilityParameters<double>>& graph, const std::string& directory) {
              int ioflags   = mio::IOF_None;
              auto ioresult = mio::write_graph<double, Model>(graph, directory, ioflags);
          });
}

} // namespace pymio

#endif //PYMIO_MOBILITY_IO_H
