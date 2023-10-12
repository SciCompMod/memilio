/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef EPI_IO_HDF5_CPP_H
#define EPI_IO_HDF5_CPP_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_HDF5

#include "memilio/utils/compiler_diagnostics.h"

MSVC_WARNING_DISABLE_PUSH(4268 4251)

#include "hdf5.h"

MSVC_WARNING_POP()

namespace mio
{

/**
 * RAII for HDF5 file handles.
 */
struct H5File {
    hid_t id;
    ~H5File()
    {
        H5Fclose(id);
    }
};

/**
 * RAII for HDF5 group handles.
 */
struct H5Group {
    hid_t id;
    ~H5Group()
    {
        H5Gclose(id);
    }
};

/**
 * RAII for HDF5 data space handles.
 */
struct H5DataSpace {
    hid_t id;
    ~H5DataSpace()
    {
        H5Sclose(id);
    }
};

/**
 * RAII for HDF5 data set handles.
 */
struct H5DataSet {
    hid_t id;
    ~H5DataSet()
    {
        H5Dclose(id);
    }
};

/**
 * Verifies a return value from the HDF5 C API.
 * Uses mio::failure to report an error if the value (the first macro argument) is negative,
 * so it must be used only in functions that return IOResult.
 * All additional macro arguments after the first are passed to an overload of mio::failure.
 * Examples:
 * MEMILIO_H5_CHECK(id, StatusCode::FileNotFound); // returns mio::failure(StatusCode::FileNotFound) if error
 * MEMILIO_H5_CHECK(id, StatusCode::FileNotFound, filename) // returns mio::failure(StatusCode::FileNotFound, filename) if error
 */
#define MEMILIO_H5_CHECK(id, ...)                                                                                      \
    do {                                                                                                               \
        if (id < 0) {                                                                                                  \
            return ::mio::failure(__VA_ARGS__);                                                                        \
        }                                                                                                              \
    } while (0)

} // namespace mio

#endif //MEMILIO_HAS_HDF5

#endif //EPI_IO_HDF5_CPP_H
