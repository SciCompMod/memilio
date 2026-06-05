/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Rene Schmieding, Julia Bicker
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
#ifndef MIO_IO_DIRECTORIES_H
#define MIO_IO_DIRECTORIES_H

#include "memilio/config.h" // IWYU pragma: keep

#include <filesystem>
#include <string>

namespace mio
{

/**
 * @brief Returns the absolute path to the project directory.
 */
std::filesystem::path base_dir();

/**
 * @brief Returns the absolute path to the project directory.
 */
[[maybe_unused]] std::filesystem::path data_dir();

/**
 * @brief Returns the absolute path to a common ouput directory for the code examples.
 * The directory lies in base_dir() and has 
 * @param[in] example_name Name of the example (e.g. the filename without .cpp).
 */
[[maybe_unused]] std::filesystem::path example_results_dir(const std::string& example_name);

} // namespace mio

#endif // MIO_IO_DIRECTORIES_H
