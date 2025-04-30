/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/timer/basic_timer.h"
#include "memilio/timer/definitions.h"
#include "memilio/timer/timers.h"
#include "gtest/gtest.h"

TEST(TimingTest, BasicTimer)
{
    // test the basic functionality of BasicTimer
    mio::timing::BasicTimer bt{};
    EXPECT_EQ(bt.get_elapsed_time(), mio::timing::DurationType{0});
    bt.start();
    bt.stop();
    // since the timers are not very precise (though accurate), we simply test that *any* measurement was taken
    EXPECT_GT(bt.get_elapsed_time(), mio::timing::DurationType{0});
}

TEST(TimingTest, AutoTimer)
{
    // check that AutoTimer correctly starts and stops a timer
    mio::timing::BasicTimer bt{};
    EXPECT_EQ(bt.get_elapsed_time(), mio::timing::DurationType{0});
    {
        mio::timing::AutoTimer at(bt);
    }
    // since the timers are not very precise (though accurate), we simply test that *any* measurement was taken
    EXPECT_GT(bt.get_elapsed_time(), mio::timing::DurationType{0});
}
