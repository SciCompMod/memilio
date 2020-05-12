#include <gtest/gtest.h>
#include <epidemiology/logging.h>

int main(int argc, char** argv)
{
    epi::set_log_level(epi::LogLevel::warn);
    ::testing::InitGoogleTest(&argc, argv);
    int retval = RUN_ALL_TESTS();
    return retval;
}
