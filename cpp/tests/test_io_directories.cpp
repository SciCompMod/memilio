/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Rene Schmieding
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
#include "memilio/io/directories.h"
#include "gtest/gtest.h"

TEST(TestDirectories, base_dir)
{
    auto base_dir = mio::base_dir();
    // check that the path exists
    EXPECT_TRUE(std::filesystem::exists(base_dir));
    EXPECT_TRUE(std::filesystem::is_directory(base_dir));
    // check that the path is correct, by sampling some fixed paths from project files
    EXPECT_TRUE(std::filesystem::exists(base_dir / "cpp" / "memilio"));
    EXPECT_TRUE(std::filesystem::exists(base_dir / "pycode" / "memilio-simulation"));
}

TEST(TestDirectories, data_dir)
{
    auto data_dir = mio::data_dir();
    // check that the path or its parent exists
    EXPECT_TRUE(std::filesystem::is_directory(data_dir) ||
                (data_dir.has_parent_path() && std::filesystem::is_directory(data_dir.parent_path())));
    // no assumptions on contents
}

TEST(TestDirectories, example_results_dir)
{
    auto ex_name  = ".__hidden_test_directory__";
    auto exrs_dir = mio::example_results_dir(ex_name);
    // check that the path does *not* exist, as examples are expected to create their own directories
    EXPECT_FALSE(std::filesystem::exists(exrs_dir));
    // check composition
    EXPECT_EQ(mio::base_dir() / "example_results" / ex_name, exrs_dir);
}
