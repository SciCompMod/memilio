/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Julia Bicker
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
#ifndef MIO_BASE_DIR_H
#define MIO_BASE_DIR_H

#include "memilio/config.h"

#include <string>

namespace mio
{

/**
 * @brief Returns the absolute path to the project directory.
 */
const static std::string base_dir()
{
    return MEMILIO_BASE_DIR;
}

} // namespace mio

#endif // MIO_BASE_DIR_H
