/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Julia Bicker, Rene Schmieding
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
#ifndef MIO_UTILS_DIRECTORIES_H
#define MIO_UTILS_DIRECTORIES_H

#include "memilio/config.h" // IWYU pragma: keep
#include "memilio/utils/stl_util.h"

#include <string>

namespace mio
{

/**
 * @brief Returns the absolute path to the project directory.
 */
const static std::string base_dir()
{
    return details::MEMILIO_BASE_DIR;
}

/**
 * @brief Returns the absolute path to the project directory.
 */
[[maybe_unused]] const static std::string data_dir()
{
    return details::MEMILIO_DATA_DIR;
}

/**
 * @brief Returns the absolute path to a common ouput directory for the code examples.
 */
[[maybe_unused]] const static std::string example_results_dir(const std::string& example_name)
{
    // the last empty string is used to end the output path in a /
    const static std::string dir = path_join(base_dir(), "example_results", example_name, "");
    return dir;
}

} // namespace mio

#endif // MIO_UTILS_DIRECTORIES_H
