#include <epidemiology/utils/logging.h>
#include <epidemiology/utils/random_number_generator.h>
#include <gtest/gtest.h>

#if HAVE_EPI_IO
#include <epidemiology_io/io.h>
#endif

int main(int argc, char** argv)
{
    epi::set_log_level(epi::LogLevel::warn);
    epi::log_thread_local_rng_seeds(epi::LogLevel::warn);

    ::testing::InitGoogleTest(&argc, argv);
    int retval = RUN_ALL_TESTS();
    return retval;
}
