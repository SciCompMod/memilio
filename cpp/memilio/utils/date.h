/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#ifndef EPI_UTILS_DATE_H
#define EPI_UTILS_DATE_H

#include "memilio/io/io.h"
#include "memilio/utils/logging.h"
#include <string>
#include <iostream>
#include <tuple>
#include <array>
#include <numeric>
#include <algorithm>
#include <cassert>

namespace mio
{

/**
 * @brief Computes if a given year is a leap year
 * @param year year as integer
 * @return true if year is a leap year, false otherwise
 */
inline bool is_leap_year(int year)
{
    return ((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0);
}

/**
 * Simple date representation as year, month, and day.
 * month in [1, 12], day in [1, 31]
 */
struct Date {
    /**
     * default constructor.
     */
    Date() = default;

    /**
     * initializing constructor.
     * @param y year as int
     * @param m month as int, 1 - 12
     * @param d day as int, 1 - 31
     */
    Date(int y, int m, int d)

        : year(y)
        , month(m)
        , day(d)
    {
        assert(month > 0 && month < 13);
        assert(day > 0 &&
               day <= ((is_leap_year(year)) ? month_lengths_leap_year[month - 1] : month_lengths[month - 1]));
    }

    /**
     * equality comparison operators.
     * @param other another date.
     */
    //@{
    bool operator==(const Date& other) const
    {
        return year == other.year && month == other.month && day == other.day;
    }
    bool operator!=(const Date& other) const
    {
        return !(*this == other);
    }
    bool operator<(const Date& other) const
    {
        if (std::tie(year, month, day) < std::tie(other.year, other.month, other.day)) {
            return true;
        }
        else {
            return false;
        }
    }
    bool operator<=(const Date& other) const
    {
        return !(other < *this);
    }
    bool operator>(const Date& other) const
    {
        return other < *this;
    }
    bool operator>=(const Date& other) const
    {
        return !(*this < other);
    }
    //@}

    /**
     * Formats the date into a string in ISO 8601 format (YYYY-MM-DD).
     * @return A string representing the date in ISO 8601 format.
     */
    std::string to_iso_string() const
    {
        // the format after ":" reads as
        // 1) '0' -> fill with zeros
        // 2) '>' -> align text right
        // 3) '4' or '2' -> specify the width (4 for year, 2 for month and day)
        return fmt::format("{:0>4}-{:0>2}-{:0>2}", year, month, day);
    }

    /**
     * Overload for stream operator to use the ISO 8601 format.
     * @param os Output stream.
     * @param date Date to output.
     * @return Reference to the output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Date& date)
    {
        return os << date.to_iso_string();
    }

    /**
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Date");
        obj.add_element("Year", year);
        obj.add_element("Month", month);
        obj.add_element("Day", day);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Date> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Date");
        auto y   = obj.expect_element("Year", Tag<int>{});
        auto m   = obj.expect_element("Month", Tag<int>{});
        auto d   = obj.expect_element("Day", Tag<int>{});
        return apply(
            io,
            [](auto&& y_, auto&& m_, auto&& d_) -> IOResult<Date> {
                if (m_ <= 0 || m_ > 12)
                    return failure(StatusCode::OutOfRange, "Month must be between 1 and 12 (inclusive).");
                if (d_ <= 0 || d_ > ((is_leap_year(y_)) ? month_lengths_leap_year[m_ - 1] : month_lengths[m_ - 1]))
                    return failure(StatusCode::OutOfRange, "Day is not valid for the given month.");
                return success(Date{y_, m_, d_});
            },
            y, m, d);
    }

    int year;
    int month;
    int day;
    static constexpr std::array<int, 12> month_lengths           = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    static constexpr std::array<int, 12> month_lengths_leap_year = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
};

/**
 * @brief Format date objects using the ISO notation for logging with spdlog.
 * @param d date object.
 */
inline std::string format_as(const mio::Date& d)
{
    return d.to_iso_string();
}

/**
 * @brief Computes the length of a month for a given date.
 * @param date date.
 * @return length of month for given date
 */
inline int get_month_length(Date date)
{
    return ((is_leap_year(date.year)) ? date.month_lengths_leap_year[date.month - 1]
                                      : date.month_lengths[date.month - 1]);
}

/**
 * @brief Computes the cumulative number of days at the end of each month for a given year.
 * @param date Date object representing the year we use to compute the cumulative days.
 * @return array with partial sum for each month
 */
inline std::array<int, 12> calculate_partial_sum_of_months(const Date& date)
{
    std::array<int, 12> part_sum;
    if (is_leap_year(date.year)) {
        std::partial_sum(date.month_lengths_leap_year.begin(), date.month_lengths_leap_year.end(), part_sum.begin());
    }
    else {
        std::partial_sum(date.month_lengths.begin(), date.month_lengths.end(), part_sum.begin());
    }
    return part_sum;
}

/**
 * parses a date from a string.
 * uses fixed format YYYY.MM.DD or YYYY-MM-DD.
 * @param date_str date as a string.
 * @return parsed date.
 */
inline IOResult<Date> parse_date(const std::string& date_str)
{
    try {
        Date date;
        date.year  = std::stoi(date_str.substr(0, 4));
        date.month = std::stoi(date_str.substr(5, 2));
        date.day   = std::stoi(date_str.substr(8, 2));
        if (date.month < 1 || date.month > 12 || date.day < 1 || date.day > get_month_length(date)) {
            return failure(StatusCode::OutOfRange, "Argument is not a valid date.");
        }
        return success(date);
    }
    catch (const std::invalid_argument&) {
        return failure(StatusCode::InvalidValue, "Argument ist not a valid date string.");
    }
}

/**
 * @brief Computes the new date corresponding to a given date and a offset in days.
 * @param date date.
 * @param offset_days offset in days.
 * @return new date that is date + offset_days.
 */
inline Date offset_date_by_days(Date date, int offset_days)
{
    auto year  = date.year;
    auto month = date.month;
    auto day   = date.day;
    assert(month > 0 && month < 13 && day > 0 && day <= get_month_length(date));

    if (day + offset_days > 0 && day + offset_days <= get_month_length(date)) {
        return {year, month, day + offset_days};
    }
    else {
        auto part_sum = calculate_partial_sum_of_months(date);

        int day_in_year = day + offset_days;
        if (month > 1) {
            // take month-2 since end of last month has to be found and due to start at 0 of C++ (against January=1 in date format)
            day_in_year += part_sum[month - 2];
        }

        if (day_in_year > 0 && day_in_year <= part_sum[11]) {
            auto iter = std::find_if(part_sum.begin(), part_sum.end(), [day_in_year](auto s) {
                return day_in_year <= s;
            });
            int i     = static_cast<int>(iter - part_sum.begin());
            return {year, i + 1, day_in_year - (i > 0 ? part_sum[static_cast<size_t>(i - 1)] : 0)};
        }
        else {
            if (day_in_year > 0) {
                return offset_date_by_days({year + 1, 1, 1}, day_in_year - part_sum[11] - 1);
            }
            else {
                return offset_date_by_days({year - 1, 12, 31}, day_in_year);
            }
        }
    }
}

/**
 * @brief Computes the day in year based on a given date.
 * @param date date
 * @return day in year, starting January, 1st, with 1.
 */
inline int get_day_in_year(Date date)
{
    auto month = date.month;
    auto day   = date.day;
    assert(month > 0 && month < 13 && day > 0 && day <= get_month_length(date));

    if (month > 1) {

        auto part_sum = calculate_partial_sum_of_months(date);

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
 * @param date1 first date.
 * @param date2 second date.
 * @return offset in days between the two dates.
 */
inline int get_offset_in_days(Date date1, Date date2)
{
    auto year1  = date1.year;
    auto month1 = date1.month;
    auto day1   = date1.day;

    auto year2  = date2.year;
    auto month2 = date2.month;
    auto day2   = date2.day;

    if (year1 == year2 && month1 == month2) {
        return day1 - day2;
    }
    else {
        int day_in_year1 = get_day_in_year(date1);
        int day_in_year2 = get_day_in_year(date2);

        if (year1 < year2) {
            int sum_days = 0;
            for (int i = year1; i < year2; i++) {
                sum_days += 365;
                if (is_leap_year(i)) {
                    sum_days += 1;
                }
            }
            return -(sum_days - day_in_year1) - day_in_year2;
        }
        else if (year1 > year2) {
            int sum_days = 0;
            for (int i = year2; i < year1; i++) {
                sum_days += 365;
                if (is_leap_year(i)) {
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

} // end namespace mio

#endif // EPI_UTILS_DATE_H
