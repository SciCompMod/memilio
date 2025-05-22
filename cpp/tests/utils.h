/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef MIO_TEST_UTILS_H
#define MIO_TEST_UTILS_H

#include "gtest/gtest.h"
#include "memilio/config.h"

namespace mio
{

/**
 * @brief Configures the mode of death tests based on the OpenMP configuration
 */
inline void set_death_test_mode()
{
#ifdef MEMILIO_ENABLE_OPENMP
    // If OpenMP is enabled, use "threadsafe" mode to silence gtest warnings.
    GTEST_FLAG_SET(death_test_style, "threadsafe");
#else
    // Without OpenMP, use "fast" mode to avoid GCov / sanitizer related errors.
    GTEST_FLAG_SET(death_test_style, "fast");
#endif
}

} // namespace mio
#endif //MIO_TEST_UTILS_H
