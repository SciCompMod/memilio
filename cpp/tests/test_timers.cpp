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
#include "utils.h"
#include "memilio/timer/auto_timer.h"
#include "memilio/timer/list_printer.h"
#include "memilio/timer/table_printer.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <utility>

struct TimingTest : public ::testing::Test {
    static void SetUpTestSuite()
    {
        // avoid uneccessary prints
        mio::timing::TimerRegistrar::get_instance().disable_final_timer_summary();
        // named timers are not reset when using the gtest_repeat flag, so we reset them manually here
        mio::timing::NamedTimer<"TestTimer", "TimingTest::AutoTimer">::get_instance().reset();
        mio::timing::NamedTimer<"TestTimer", "TimingTest::NamedTimer">::get_instance().reset();
    }

private:
    inline static mio::timing::BasicTimer m_timer;

public:
    static std::list<mio::timing::TimerRegistration> printable_timers()
    {
        return {{"A", "A", m_timer, 0},
                {"B", "", m_timer, 0},
                {"B", "", m_timer, 1},
                {"C has a fairly long name", "C", m_timer, 0},
                {"D", "D", m_timer, 2}};
    }

    // something to make sure the timers do not measure 0
    static void wait_briefly()
    {
        std::this_thread::sleep_for(std::chrono::nanoseconds(100));
    }
};

// very simple printer that only outputs the number of registrations
struct SizePrinter : public mio::timing::Printer {
    void print(const std::list<mio::timing::TimerRegistration>& list, std::ostream& out) override
    {
        out << list.size();
    }
};

TEST_F(TimingTest, BasicTimer)
{
    // replace spdlog::default_logger
    mio::RedirectLogger logger;
    logger.capture();
    // test the basic functionality of BasicTimer
    mio::timing::BasicTimer bt{};
    EXPECT_EQ(bt.get_elapsed_time(), mio::timing::DurationType{0});
    bt.start();
    wait_briefly();
    bt.stop();
    // since the timers are not very precise (though accurate), we simply test that *any* measurement was taken
    EXPECT_GT(bt.get_elapsed_time(), mio::timing::DurationType{0});
    bt.reset();
    // check that reset works
    EXPECT_EQ(bt.get_elapsed_time(), mio::timing::DurationType{0});
    // check that there are no error reports
    EXPECT_EQ(logger.read(), "");
#ifndef NDEBUG
    // almost all of BasicTimer, including its destructor, expects to be called while stopped - except for stop itself
    bt.start();
    bt.start();
    EXPECT_THAT(logger.read(), ::testing::HasSubstr("expected to be stopped."));
    bt.get_elapsed_time();
    EXPECT_THAT(logger.read(), ::testing::HasSubstr("expected to be stopped."));
    bt.reset();
    EXPECT_THAT(logger.read(), ::testing::HasSubstr("expected to be stopped."));
    bt.stop();
    bt.stop();
    EXPECT_THAT(logger.read(), ::testing::HasSubstr("expected to be started."));
    { // trigger destructor
        mio::timing::BasicTimer{}.start();
    }
    EXPECT_THAT(logger.read(), ::testing::HasSubstr("expected to be stopped."));
#endif
    // restore spdlog::default_logger
    logger.release();
}

TEST_F(TimingTest, NamedTimer)
{
    // check that getters for instance, name and scope work
    // Note: named timers have to be reset in SetUpTestSuite
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
    // Note: named timers have to be reset in SetUpTestSuite
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
    auto& reg                = mio::timing::TimerRegistrar::get_instance();
    const int old_num_timers = (int)reg.get_register().size(); // this may be bigger than 0, depending on other tests
    std::vector<mio::timing::BasicTimer> bts(mio::omp::get_max_threads());

    // check add_timer by registering a timer (per omp thread)
    PRAGMA_OMP(parallel num_threads(mio::omp::get_max_threads()))
    {
        const auto i = mio::omp::get_thread_id();
        reg.add_timer({"name", "", bts[i], i}); // reusing names is bad. but the registrar does not care
    }
    // -> verify by register size
    const int new_num_timers = old_num_timers + mio::omp::get_max_threads();
    ASSERT_EQ(reg.get_register().size(), new_num_timers);
    // -> verify by registered thread_id (the last registered timers should be 0,...,num_threads, in any order)
    std::vector<bool> ids_found(mio::omp::get_max_threads(), false);
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

TEST_F(TimingTest, ListPrinter)
{
    // test printer against a known good configuration. changes to the printer will break this test. in that case,
    // manually verify the new behaviour, and update the reference_output.

    const std::string reference_output = "All Timers: 5\n"
                                         "  A::A: 0.000000e+00 (0)\n"
                                         "  B: 0.000000e+00 (0)\n"
                                         "  B: 0.000000e+00 (1)\n"
                                         "  C::C has a fairly long name: 0.000000e+00 (0)\n"
                                         "  D::D: 0.000000e+00 (2)\n"
                                         "Unique Timers (accumulated): 4\n"
                                         "  A::A: 0.000000e+00\n"
                                         "  B: 0.000000e+00\n"
                                         "  C::C has a fairly long name: 0.000000e+00\n"
                                         "  D::D: 0.000000e+00\n";
    std::stringstream out;

    mio::timing::ListPrinter p;
    p.print(this->printable_timers(), out);

    EXPECT_EQ(out.str(), reference_output);
}

TEST_F(TimingTest, TablePrinter)
{
    // test printer against a known good configuration. changes to the printer will break this test. in that case,
    // manually verify the new behaviour, and update the reference_output.

    {
        const std::string reference_output =
            "+-----------------------------+--------------+--------------+--------------+--------------+----------+\n"
            "| Timers                      | Elapsed Time | Min          | Max          | Average      | #Threads |\n"
            "+-----------------------------+--------------+--------------+--------------+--------------+----------+\n"
            "| A::A                        | 0.000000e+00 |           -- |           -- |           -- |        1 |\n"
            "| B                           | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 |        2 |\n"
            "| C::C has a fairly long name | 0.000000e+00 |           -- |           -- |           -- |        1 |\n"
            "| D::D                        | 0.000000e+00 |           -- |           -- |           -- |        1 |\n"
            "+-----------------------------+--------------+--------------+--------------+--------------+----------+\n";

        std::stringstream out;

        mio::timing::TablePrinter p;
        p.print(this->printable_timers(), out);

        EXPECT_EQ(out.str(), reference_output);
    }

    // same as above, but change the format
    {
        const std::string reference_output =
            "+-----------------------------+--------------+--------+--------+---------+----------+\n"
            "| Timers                      | Elapsed Time | Min    | Max    | Average | #Threads |\n"
            "+-----------------------------+--------------+--------+--------+---------+----------+\n"
            "| A::A                        |         0.00 |     -- |     -- |      -- |        1 |\n"
            "| B                           |         0.00 |   0.00 |   0.00 |    0.00 |        2 |\n"
            "| C::C has a fairly long name |         0.00 |     -- |     -- |      -- |        1 |\n"
            "| D::D                        |         0.00 |     -- |     -- |      -- |        1 |\n"
            "+-----------------------------+--------------+--------+--------+---------+----------+\n";

        std::stringstream out;

        mio::timing::TablePrinter p;
        p.set_time_format("{:6.2f}");
        p.print(this->printable_timers(), out);

        EXPECT_EQ(out.str(), reference_output);
    }
}
