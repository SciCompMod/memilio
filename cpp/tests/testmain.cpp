/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Martin Siggel
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
#include "memilio/utils/logging.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include <gtest/gtest.h>

int main(int argc, char** argv)
{
    mio::mpi::init();
    mio::set_log_level(mio::LogLevel::info); // only used for logging testing environment info
    mio::log_thread_local_rng_seeds(mio::LogLevel::info);

#ifndef NDEBUG
    // in debug builds:
    // suppress output from successful tests, so output from failures is easy to find
    GTEST_FLAG_SET(brief, true);
    // shuffle order of tests, to make sure they are not interdependant
    GTEST_FLAG_SET(shuffle, true);
#endif
    ::testing::InitGoogleTest(&argc, argv);

    mio::set_log_level(mio::LogLevel::warn); // main log level for testing
    int retval = RUN_ALL_TESTS();
    mio::mpi::finalize();
    return retval;
}
