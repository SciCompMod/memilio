/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Wadim Koslow
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
#ifndef IO_H
#define IO_H

#include "epidemiology/utils/io.h"
#include <string>

namespace epi
{

/**
 * @brief Returns the current working directory name
 */ 
std::string get_current_dir_name();

/**
 * @brief Creates a directory on the file system
 * returns true if the creation was successful.
 */
IOResult<bool> create_directory(std::string const& rel_path, std::string& abs_path);

bool file_exists(std::string const& rel_path, std::string& abs_path);

} // namespace epi

#endif // IO_H
