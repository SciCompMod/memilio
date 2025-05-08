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
#include "memilio/timer/auto_timer.h"
#include "memilio/utils/mioomp.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

struct TimingTest : public ::testing::Test {
    static void SetUpTestSuite()
    {
        mio::timing::TimerRegistrar::get_instance().disable_final_timer_summary();
    }
};

struct SizePrinter : public mio::timing::Printer {
    void print(const std::list<mio::timing::TimerRegistration>& list, std::ostream& out) override
    {
        out << list.size();
    }
};

void wait_briefly()
{
    std::this_thread::sleep_for(std::chrono::nanoseconds(100));
}

TEST_F(TimingTest, BasicTimer)
{
    // test the basic functionality of BasicTimer
    mio::timing::BasicTimer bt{};
    EXPECT_EQ(bt.get_elapsed_time(), mio::timing::DurationType{0});
    bt.start();
    wait_briefly();
    bt.stop();
    // since the timers are not very precise (though accurate), we simply test that *any* measurement was taken
    EXPECT_GT(bt.get_elapsed_time(), mio::timing::DurationType{0});
}

TEST_F(TimingTest, NamedTimer)
{
    // check that getters for instance, name and scope work
    auto& nt = mio::timing::NamedTimer<"TestTimer", "TimingTest::NamedTimer">::get_instance();
    EXPECT_STREQ("TestTimer", nt.name().c_str());
    EXPECT_STREQ("TimingTest::NamedTimer", nt.scope().c_str());
    // timing functionality is inherited from BasicTimer, but we check it again anyways
    EXPECT_EQ(nt.get_elapsed_time(), mio::timing::DurationType{0});
    nt.start();
    wait_briefly();
    nt.stop();
    // since the timers are not very precise (though accurate), we simply test that *any* measurement was taken
    EXPECT_GT(nt.get_elapsed_time(), mio::timing::DurationType{0});
}

TEST_F(TimingTest, AutoTimer)
{
    // check that AutoTimer correctly starts and stops a timer by reference
    mio::timing::BasicTimer bt{};
    EXPECT_EQ(bt.get_elapsed_time(), mio::timing::DurationType{0});
    {
        mio::timing::AutoTimer at(bt);
        wait_briefly();
    }
    auto& nt = mio::timing::NamedTimer<"TestTimer", "TimingTest::AutoTimer">::get_instance();
    EXPECT_EQ(nt.get_elapsed_time(), mio::timing::DurationType{0});
    {
        mio::timing::AutoTimer<"TestTimer", "TimingTest::AutoTimer"> at;
        wait_briefly();
    }
    // since the timers are not very precise (though accurate), we simply test that *any* measurement was taken
    EXPECT_GT(bt.get_elapsed_time(), mio::timing::DurationType{0});
    EXPECT_GT(nt.get_elapsed_time(), mio::timing::DurationType{0});
}

TEST_F(TimingTest, TimerRegistrar)
{
    // check member functions by using them for adding basic timers:
    // use get_instance and get_register. implicitly verified by later tests
    auto& reg                 = mio::timing::TimerRegistrar::get_instance();
    const auto old_num_timers = reg.get_register().size(); // this may be bigger than 0, depending on other tests
    std::vector<mio::timing::BasicTimer> bts(mio::get_omp_num_threads());

    // check add_timer by registering a timer (per omp thread)
    PRAGMA_OMP(parallel)
    {
        const auto i = mio::get_omp_thread_id();
        reg.add_timer({"name", "", bts[i], i}); // reusing names is bad. but the registrar does not care
    }
    // -> verify by register size
    const auto new_num_timers = old_num_timers + mio::get_omp_num_threads();
    ASSERT_EQ(reg.get_register().size(), new_num_timers);
    // -> verify by registered thread_id (the last registered timers should be 0,...,num_threads, in any order)
    std::vector<bool> ids_found(mio::get_omp_num_threads(), false);
    auto reg_itr = reg.get_register().crbegin();
    for (size_t i = 0; i < ids_found.size(); i++) {
        ids_found[reg_itr->thread_id] = true;
        ++reg_itr;
    }
    EXPECT_TRUE(std::all_of(ids_found.begin(), ids_found.end(), [](auto&& b) {
        return b;
    }));

    // check set_printer. verification is done together with print_timers
    auto p = std::make_unique<SizePrinter>();
    reg.set_printer(std::move(p));
    // check print_timers by using the SizePrinter (defined above), which only writes reg.get_register().size()
    int size = -1;
    std::stringstream out;
    reg.print_timers(out);
    out >> size;
    EXPECT_EQ(size, new_num_timers);
}

TEST_F(TimingTest, qualified_name)
{
    // test that name and scope are concatonated as expected
    EXPECT_EQ(mio::timing::qualified_name("name", ""), "name");
    EXPECT_EQ(mio::timing::qualified_name("name", "scope"), "scope::name");
}
