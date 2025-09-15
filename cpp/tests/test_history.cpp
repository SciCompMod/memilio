/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding, Sascha Korf
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
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "memilio/io/history.h"

class example
{
public:
    int a            = 1;
    int b            = 2;
    int current_time = 0;
};

struct LogPair : mio::LogAlways {
    using Type = std::pair<int, int>;
    static Type log(const example& ex)
    {
        return {ex.a, ex.b};
    }
};
struct LogAOnce : mio::LogOnce {
    using Type = int;
    static Type log(const example& ex)
    {
        return ex.a;
    }
};

struct LogStepIf {
    using Type = int;
    static Type log(const example& ex)
    {
        return ex.current_time;
    }
    static bool should_log(const example& ex)
    {
        return ex.current_time == 0;
    }
};

TEST(HistoryObject, log)
{
    example ex;
    mio::HistoryWithMemoryWriter<LogPair, LogAOnce, LogStepIf> history;
    int n_runs = 2;
    for (int i = 0; i < n_runs; i++) {
        ex.current_time = i;
        history.log(ex);
    }
    auto data = history.get_log();
    ASSERT_EQ(std::get<0>(data).size(), n_runs);
    ASSERT_EQ(std::get<1>(data).size(), 1); //cause we log this only once
    ASSERT_EQ(std::get<0>(data)[0], std::get<0>(data)[1]);
    ASSERT_EQ(std::get<0>(data)[0], std::make_pair(ex.a, ex.b));
    ASSERT_EQ(std::get<1>(data)[0], ex.a);
    ASSERT_EQ(std::get<2>(data)[0], 0);
}
