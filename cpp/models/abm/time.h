/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef EPI_ABM_TIME_H
#define EPI_ABM_TIME_H

namespace mio
{

/**
 * a duration of time.
 * resolution 1 second.
 */
class TimeSpan
{
public:
    /**
     * default ctor, uninitialized.
     */
    TimeSpan() = default;

    /**
     * creates a TimeSpan that represents a number of seconds.
     */
    explicit TimeSpan(int seconds)
        : m_seconds(seconds)
    {
    }

    /**
     * length of time in days.
     */
    double days() const
    {
        return double(m_seconds) / (24 * 60 * 60);
    }

    /**
     * length of time in hours.
     */
    double hours() const
    {
        return double(m_seconds) / (60 * 60);
    };
    
    /**
     * length of time in seconds.
     */
    int seconds() const
    {
        return m_seconds;
    }

    /**
     * comparison operators
     * @{
     */
    bool operator==(const TimeSpan& other) const
    {
        return m_seconds == other.m_seconds;
    }
    bool operator!=(const TimeSpan& other) const
    {
        return !(*this == other);
    }
    bool operator<(const TimeSpan& other) const
    {
        return m_seconds < other.m_seconds;
    }
    bool operator<=(const TimeSpan& other) const
    {
        return m_seconds <= other.m_seconds;
    }
    bool operator>(const TimeSpan& other) const
    {
        return m_seconds > other.m_seconds;
    }
    bool operator>=(const TimeSpan& other) const
    {
        return m_seconds >= other.m_seconds;
    }
    /**@}*/

    /**
     * numeric operators for addition, subtraction, and scalar integer multiplication and division. 
     * @{
     */
    TimeSpan operator+(const TimeSpan& s) const
    {
        return TimeSpan{m_seconds + s.m_seconds};
    }
    TimeSpan& operator+=(const TimeSpan& s)
    {
        m_seconds += s.m_seconds;
        return *this;
    }
    TimeSpan operator-(const TimeSpan& s) const
    {
        return TimeSpan{m_seconds - s.m_seconds};
    }
    TimeSpan& operator-=(const TimeSpan& s)
    {
        m_seconds -= s.m_seconds;
        return *this;
    }

    TimeSpan operator*(int f) const
    {
        return TimeSpan{m_seconds * f};
    }
    TimeSpan& operator*=(int f)
    {
        m_seconds *= f;
        return *this;
    }
    TimeSpan operator/(int f) const
    {
        return TimeSpan{m_seconds / f};
    }
    TimeSpan& operator/=(int f)
    {
        m_seconds /= f;
        return *this;
    }
    /**@}*/

private:
    int m_seconds;
};

/**
 * represents a point in time.
 * seconds from an unspecified monday at 00:00 (epoch). 
 */
class TimePoint
{
public:
    /**
     * default ctor, unitinialized.
     */
    TimePoint() = default;
    /**
     * creates a TimePoint from a specified number of seconds.
     */
    explicit TimePoint(int seconds)
        : m_seconds(seconds)
    {
    }

    /**
     * time since the epoch in days.
     */
    double days() const
    {
        return double(m_seconds) / (24 * 60 * 60);
    }
    /**
     * time since the epoch in hours.
     */
    double hours() const
    {
        return double(m_seconds) / (60 * 60);
    };    
    /**
     * time since the epoch in seconds.
     */
    int seconds() const
    {
        return m_seconds;
    }

    /**
     * index of current day of the week (0,...,6 = Mo,...,Sun)
     */
    int day_of_week() const
    {
        return int(days()) % 7;
    }

    /**
     * hour in the current day (0 - 23).
     */
    int hour_of_day() const
    {
        return int(hours()) % 24;
    }

    /**
     * time since midnight
     */
    TimeSpan time_since_midnight() const
    {
        return  TimeSpan(seconds() - ((int)days())*60*60*24);
    }

    /**
     * comparison operators.
     * @{
     */
    bool operator==(const TimePoint& other) const
    {
        return m_seconds == other.m_seconds;
    }
    bool operator!=(const TimePoint& other) const
    {
        return !(*this == other);
    }
    bool operator<(const TimePoint& other) const
    {
        return m_seconds < other.m_seconds;
    }
    bool operator<=(const TimePoint& other) const
    {
        return m_seconds <= other.m_seconds;
    }
    bool operator>(const TimePoint& other) const
    {
        return m_seconds > other.m_seconds;
    }
    bool operator>=(const TimePoint& other) const
    {
        return m_seconds >= other.m_seconds;
    }
    /**@}*/

    /**
     * add or subtract a TimeSpan.
     * @{
     */
    TimePoint operator+(const TimeSpan& s) const
    {
        return TimePoint{m_seconds + s.seconds()};
    }
    TimePoint& operator+=(const TimeSpan& s)
    {
        m_seconds += s.seconds();
        return *this;
    }
    TimePoint operator-(const TimeSpan& s) const
    {
        return TimePoint{m_seconds - s.seconds()};
    }
    TimePoint& operator-=(const TimeSpan& s)
    {
        m_seconds -= s.seconds();
        return *this;
    }
    /**@}*/

    /**
     * TimeSpan difference between two time points.
     */
    TimeSpan operator-(const TimePoint& p2)
    {
        return TimeSpan{m_seconds - p2.seconds()};
    }

private:
    int m_seconds;
};

/**
 * create a TimeSpan of a specified number of seconds.
 * @param seconds number of seconds in the time span.
 */
inline TimeSpan seconds(int seconds)
{
    return TimeSpan(seconds);
}

/**
 * create a TimeSpan of a specified number of minutes.
 * @param minutes number of minutes in the time span.
 */
inline TimeSpan minutes(int minutes)
{
    return TimeSpan(minutes * 60);
}


/**
 * create a TimeSpan of a specified number of hours.
 * @param seconds number of hours in the time span.
 */
inline TimeSpan hours(int hours)
{
    return TimeSpan(hours * 60 * 60);
}

/**
 * create a TimeSpan with a specified number of days.
 * @param seconds number of days in the time span.
 */
inline TimeSpan days(int days)
{
    return TimeSpan(days * 24 * 60 * 60);
}

} // namespace mio

#endif