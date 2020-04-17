#include <gtest/gtest.h>

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int retval = RUN_ALL_TESTS();
    return retval;
}
