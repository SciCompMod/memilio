/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kuehn
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
#include "memilio/io/result_io.h"

#ifdef MEMILIO_HAS_HDF5

#include "memilio/io/hdf5_cpp.h"
#include "memilio/math/eigen_util.h"
#include "memilio/epidemiology/damping.h"

#include <vector>
#include <iostream>
#include <string>

namespace mio
{

herr_t store_group_name(hid_t /*id*/, const char* name, const H5L_info_t* /*linfo*/, void* opdata)
{
    auto h5group_names = reinterpret_cast<std::vector<std::string>*>(opdata);
    h5group_names->push_back(name);
    return 0;
}

} // namespace mio

#endif //MEMILIO_HAS_HDF5
