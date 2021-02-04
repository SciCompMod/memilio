#ifndef DATE_HELPER_H
#define DATE_HELPER_H

#include <string>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <cassert>

namespace epi
{

/**
 * @brief Data format storing year, month, and day as separated int
 * @param year year of input date.
 * @param month month of input date.
 * @param day day of input date.
 * @return struct with year, month, and day of new date.
 */
struct date_sep {
    int year;
    int month;
    int day;
};

inline struct date_sep parse_date(const std::string& date_str)
{
    date_sep date;
    date.year  = std::stoi(date_str.substr(0, 4));
    date.month = std::stoi(date_str.substr(5, 2));
    date.day   = std::stoi(date_str.substr(8, 2));
    return date;
}

/**
 * @brief Computes the new date corresponding to a given date and a offset in days.
 * @param year year of input date.
 * @param month month of input date.
 * @param day day of input date.
 * @param offset offset in days.
 * @return struct with year, month, and day of new date.
 */
inline struct date_sep calc_date_offset(int year, int month, int day, int offset)
{
    assert(month > 0 && month < 13 && day > 0 && day < 32);

    std::vector<int> month_len;
    if (((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0)) {
        // leap year
        month_len = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    }
    else {
        month_len = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    }
    if (day + offset > 0 && day + offset <= month_len[month - 1]) {
        return {year, month, day + offset};
    }
    else {
        std::vector<int> part_sum(12, 0);
        std::partial_sum(month_len.begin(), month_len.end(), part_sum.begin());

        int day_in_year = day + offset;
        if (month > 1) {
            // take month-2 since end of last month has to be found and due to start at 0 of C++ (against January=1 in date format)
            day_in_year += part_sum[month - 2];
        }

        if (day_in_year > 0 && day_in_year <= part_sum[11]) {

            int i = 0;
            while (day_in_year > part_sum[i]) {
                i++;
            }
            return {year, i + 1, day_in_year - part_sum[i - 1]};
        }
        else {
            if(day_in_year > 0)
            {
                return calc_date_offset(year + 1, 1, 1, day_in_year - part_sum[11] - 1);
            }else{
                return calc_date_offset(year - 1, 12, 31, day_in_year);
            }
        }
    }
} // namespace epi

/**
 * @brief Computes the day in year based on a given date.
 * @param year year of input date.
 * @param month month of input date.
 * @param day day of input date.
 * @return day in year, starting January, 1st, with 1.
 */
inline int calc_day_in_year(int year, int month, int day)
{
    assert(month > 0 && month < 13 && day > 0 && day < 32);

    if (month > 1) {
        std::vector<int> month_len;
        if (((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0)) {
            // leap year
            month_len = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        }
        else {
            month_len = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
        }

        std::vector<int> part_sum(12, 0);
        std::partial_sum(month_len.begin(), month_len.end(), part_sum.begin());

        // take month-2 since end of last month has to be found and due to start at 0 of C++ (against January=1 in date format)
        int day_in_year = part_sum[month - 2] + day;

        return day_in_year;
    }
    else {
        return day;
    }
}

/**
 * @brief Computes the offset in days given two dates: first date minus second date.
 * @param year1 year of first date.
 * @param month1 month of first date.
 * @param day1 day of first date.
 * @param year2 year of second date.
 * @param month2 month of second date.
 * @param day2 day of second date.
 * @return offset offset in days between the two dates.
 */
inline int calc_offset_dates(int year1, int month1, int day1, int year2, int month2, int day2)
{
    if (year1 == year2 && month1 == month2) {
        return day1 - day2;
    }
    else {
        int day_in_year1 = calc_day_in_year(year1, month1, day1);
        int day_in_year2 = calc_day_in_year(year2, month2, day2);

        if (year1 < year2) {
            int sum_days = 0;
            for (int i = year1; i < year2; i++) {
                sum_days += 365;
                if (((i % 4 == 0) && (i % 100 != 0)) || (i % 400 == 0)) {
                    sum_days += 1;
                }
            }
            return -(sum_days - day_in_year1) - day_in_year2;
        }
        else if (year1 > year2) {
            int sum_days = 0;
            for (int i = year2; i < year1; i++) {
                sum_days += 365;
                if (((i % 4 == 0) && (i % 100 != 0)) || (i % 400 == 0)) {
                    sum_days += 1;
                }
            }
            return day_in_year1 + sum_days - day_in_year2;
        }
        else {

            return day_in_year1 - day_in_year2;
        }
    }
}

/* 
 * @brief checks if a given date corresponds to another date plus a given offset in days.
 * @param year1 year of first date.
 * @param month1 month of first date.
 * @param day1 day of first date.
 * @return true if the dates correspond according to the offset.
 */
inline bool check_date(int year, int month, int day, int offset, date_sep date)
{
    struct date_sep date_offset = calc_date_offset(year, month, day, offset);
    bool correct_date = date_offset.year == date.year && date_offset.month == date.month && date_offset.day == date.day;
    return correct_date;
}

/* 
 * @brief checks if a given date corresponds to another date plus a given offset in days
 * @param year1 year of first date
 * @param month1 month of first date
 * @param day1 day of first date
 * @return true if the dates correspond according to the offset
 */
inline bool check_date(int year, int month, int day, int offset, const std::string& date)
{
    return check_date(year, month, day, offset, parse_date(date));
}

} // end namespace epi

#endif // DATE_HELPER_H
