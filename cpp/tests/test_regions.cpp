#include "epidemiology/utils/regions.h"
#include "gtest/gtest.h"

TEST(TestRegions, get_holidays)
{
    auto a = epi::regions::de::get_holidays(epi::regions::de::StateId(9), epi::Date(2020, 10, 15), epi::Date(2020, 11, 15));
    ASSERT_EQ(a.size(), 1);
    ASSERT_EQ(a[0], std::make_pair(epi::Date(2020, 10, 31), epi::Date(2020, 11, 7)));

    auto b = epi::regions::de::get_holidays(epi::regions::de::StateId(3), epi::Date(2020, 7, 30), epi::Date(2020, 12, 31));
    ASSERT_EQ(b.size(), 3);
    ASSERT_EQ(b[0], std::make_pair(epi::Date(2020, 7, 16), epi::Date(2020, 8, 27)));
    ASSERT_EQ(b[1], std::make_pair(epi::Date(2020, 10, 12), epi::Date(2020, 10, 24)));
    ASSERT_EQ(b[2], std::make_pair(epi::Date(2020, 12, 23), epi::Date(2021, 1, 9)));
}

TEST(TestRegions, get_state_id)
{
    ASSERT_EQ(epi::regions::de::get_state_id(epi::regions::de::CountyId(1001)), epi::regions::de::StateId(1));
    ASSERT_EQ(epi::regions::de::get_state_id(epi::regions::de::CountyId(2000)), epi::regions::de::StateId(2));
    ASSERT_EQ(epi::regions::de::get_state_id(epi::regions::de::CountyId(5970)), epi::regions::de::StateId(5));
    ASSERT_EQ(epi::regions::de::get_state_id(epi::regions::de::CountyId(9161)), epi::regions::de::StateId(9));
}