/* 
* Copyright (C) 2020-2025 MEmilio
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

#include "abm/time.h"
#include "memilio/io/default_serialize.h"

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
    Generic,
    Antigen,
    PCR,

    Count
};

/**
* @brief The TestResult of a Person.
*/
struct TestResult {
    TimePoint time_of_testing{std::numeric_limits<int>::min()}; ///< The TimePoint when the Person performs the test.
    bool result{false}; ///< The test result.

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("TestResult").add("time_of_testing", time_of_testing).add("result", result);
    }
};

} // namespace abm
} // namespace mio

#endif
