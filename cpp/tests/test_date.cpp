#include <epidemiology/utils/date.h>
#include <gtest/gtest.h>

TEST(TestDate, calcOffset)
{
    auto new_date = epi::calc_date_offset(2020, 8, 30, 3);
    EXPECT_EQ(new_date.year, 2020);
    EXPECT_EQ(new_date.month, 9);
    EXPECT_EQ(new_date.day, 2);

    new_date = epi::calc_date_offset(2024, 8, 24, 6);
    EXPECT_EQ(new_date.year, 2024);
    EXPECT_EQ(new_date.month, 8);
    EXPECT_EQ(new_date.day, 30);

    new_date = epi::calc_date_offset(2024, 8, 24, -6);
    EXPECT_EQ(new_date.year, 2024);
    EXPECT_EQ(new_date.month, 8);
    EXPECT_EQ(new_date.day, 18);    

    new_date = epi::calc_date_offset(2020, 2, 28, 1);
    EXPECT_EQ(new_date.year, 2020);
    EXPECT_EQ(new_date.month, 2);
    EXPECT_EQ(new_date.day, 29);

    new_date = epi::calc_date_offset(2021, 2, 28, 1);
    EXPECT_EQ(new_date.year, 2021);
    EXPECT_EQ(new_date.month, 3);
    EXPECT_EQ(new_date.day, 1);

    new_date = epi::calc_date_offset(2021, 2, 3, -5);
    EXPECT_EQ(new_date.year, 2021);
    EXPECT_EQ(new_date.month, 1);
    EXPECT_EQ(new_date.day, 29);    

    new_date = epi::calc_date_offset(2021, 2, 28, 0);
    EXPECT_EQ(new_date.year, 2021);
    EXPECT_EQ(new_date.month, 2);
    EXPECT_EQ(new_date.day, 28);

    new_date = epi::calc_date_offset(2020, 12, 31, 1);
    EXPECT_EQ(new_date.year, 2021);
    EXPECT_EQ(new_date.month, 1);
    EXPECT_EQ(new_date.day, 1);  

    new_date = epi::calc_date_offset(2021, 12, 31, 1);
    EXPECT_EQ(new_date.year, 2022);
    EXPECT_EQ(new_date.month, 1);
    EXPECT_EQ(new_date.day, 1);   

    new_date = epi::calc_date_offset(2019, 12, 31, 367);
    EXPECT_EQ(new_date.year, 2021);
    EXPECT_EQ(new_date.month, 1);
    EXPECT_EQ(new_date.day, 1);  

    new_date = epi::calc_date_offset(2021, 12, 31, 366);
    EXPECT_EQ(new_date.year, 2023);
    EXPECT_EQ(new_date.month, 1);
    EXPECT_EQ(new_date.day, 1);            

    new_date = epi::calc_date_offset(2021, 1, 1, -1);
    EXPECT_EQ(new_date.year, 2020);
    EXPECT_EQ(new_date.month, 12);
    EXPECT_EQ(new_date.day, 31);  

    new_date = epi::calc_date_offset(2022, 1, 1, -1);
    EXPECT_EQ(new_date.year, 2021);
    EXPECT_EQ(new_date.month, 12);
    EXPECT_EQ(new_date.day, 31);        

    new_date = epi::calc_date_offset(2021, 1, 1, -367);
    EXPECT_EQ(new_date.year, 2019);
    EXPECT_EQ(new_date.month, 12);
    EXPECT_EQ(new_date.day, 31);  

    new_date = epi::calc_date_offset(2022, 1, 1, -366);
    EXPECT_EQ(new_date.year, 2020);
    EXPECT_EQ(new_date.month, 12);
    EXPECT_EQ(new_date.day, 31);             
}

TEST(TestDate, checkDate)
{
    auto checked_date = epi::check_date(2020, 8, 30, 3, "2020.09.02");
    EXPECT_EQ(checked_date, true);

    checked_date = epi::check_date(2021, 8, 24, 6, "2021.08.30");
    EXPECT_EQ(checked_date, true);

    checked_date = epi::check_date(2020, 2, 28, 1, "2020.02.29");
    EXPECT_EQ(checked_date, true);

    checked_date = epi::check_date(2021, 2, 28, 1, "2021.03.01");
    EXPECT_EQ(checked_date, true);

    checked_date = epi::check_date(2021, 2, 28, 0, "2021.02.28");
    EXPECT_EQ(checked_date, true);

    checked_date = epi::check_date(2020, 8, 30, 4, "2020.09.02");
    EXPECT_EQ(checked_date, false);

    checked_date = epi::check_date(2021, 8, 24, 7, "2020.08.30");
    EXPECT_EQ(checked_date, false);

    checked_date = epi::check_date(2020, 2, 28, 2, "2020.02.29");
    EXPECT_EQ(checked_date, false);

    checked_date = epi::check_date(2021, 2, 28, 0, "2021.03.01");
    EXPECT_EQ(checked_date, false);
}

TEST(TestDate, calcDayInYear)
{
    auto day = epi::calc_day_in_year(2020, 1, 21);
    EXPECT_EQ(day, 21);

    day = epi::calc_day_in_year(2020, 2, 14);
    EXPECT_EQ(day, 45);

    day = epi::calc_day_in_year(2020, 3, 2);
    EXPECT_EQ(day, 62);

    day = epi::calc_day_in_year(2020, 12, 27);
    EXPECT_EQ(day, 362);

    day = epi::calc_day_in_year(2021, 1, 21);
    EXPECT_EQ(day, 21);

    day = epi::calc_day_in_year(2021, 2, 14);
    EXPECT_EQ(day, 45);

    day = epi::calc_day_in_year(2021, 3, 2);
    EXPECT_EQ(day, 61);

    day = epi::calc_day_in_year(2021, 12, 27);
    EXPECT_EQ(day, 361);
}

TEST(TestDate, calcOffsetDates)
{
    auto offset = epi::calc_offset_dates(2020, 8, 30, 2020, 8, 15);
    EXPECT_EQ(offset, 15);

    offset = epi::calc_offset_dates(2020, 8, 30, 2020, 8, 31);
    EXPECT_EQ(offset, -1);

    offset = epi::calc_offset_dates(2020, 9, 12, 2020, 8, 30);
    EXPECT_EQ(offset, 13);

    offset = epi::calc_offset_dates(2020, 7, 25, 2020, 5, 25);
    EXPECT_EQ(offset, 61);

    offset = epi::calc_offset_dates(2021, 1, 3, 2020, 12, 31);
    EXPECT_EQ(offset, 3);

    offset = epi::calc_offset_dates(2021, 3, 3, 2020, 12, 29);
    EXPECT_EQ(offset, 64);

    offset = epi::calc_offset_dates(2021, 11, 30, 2020, 11, 30);
    EXPECT_EQ(offset, 365);

    offset = epi::calc_offset_dates(2025, 11, 30, 2020, 11, 30);
    EXPECT_EQ(offset, 5 * 365 + 1);

    offset = epi::calc_offset_dates(2019, 11, 30, 2020, 11, 30);
    EXPECT_EQ(offset, -366);
}