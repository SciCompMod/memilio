/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Khoa Nguyen
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

#ifndef MIO_ABM_TEST_TYPE_H
#define MIO_ABM_TEST_TYPE_H

#include <cstdint>

namespace mio
{
namespace abm
{

/**
 * @brief Type of a Test.
 */
enum class TestType : std::uint32_t
{
    Generic = 0,
    Antigen,
    PCR,

    Count
};

} // namespace abm
} // namespace mio

#endif