/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker, Maximilian Betz
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
#ifndef PYMIO_IO_RESULT_IO_H
#define PYMIO_IO_RESULT_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_HDF5

#include "memilio/io/result_io.h"

#include "pybind11/pybind11.h"

namespace pymio
{

template <class Model>
void bind_save_results(pybind11::module_& m)
{
    m.def("save_results",
          [&](const std::vector<std::vector<mio::TimeSeries<double>>>& ensemble_results,
              const std::vector<std::vector<Model>>& ensemble_params, const std::vector<int>& county_ids,
              const std::string& result_dir, bool save_single_runs, bool save_percentiles) {
              boost::filesystem::path dir(result_dir);
              auto ioresult = mio::save_results<double, Model>(ensemble_results, ensemble_params, county_ids, dir,
                                                               save_single_runs, save_percentiles);
              return NULL;
          });
}

} // namespace pymio

#endif // MEMILIO_HAS_HDF5

#endif //PYMIO_IO_RESULT_IO_H
