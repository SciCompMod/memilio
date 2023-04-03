/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Sascha Korf, Rene Schmieding 
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
#include <vector>
#include <tuple>
#include <iostream>
#include "memilio/utils/index.h"

//TODO: Maybe rename Log... something like Reader or something similar.
class example
{
public:
    int a = 1;
    int b = 2;
};
struct LogOnce {
};

struct LogA : LogOnce {
    using Type = int;
    static Type log(const example& ex)
    {
        return ex.a;
    }
};

struct LogAlways {
};

struct LogB : LogAlways {
    using Type = int;
    static Type log(const example& ex)
    {
        return ex.b;
    }
};

template <class... Loggers>
struct Writer {
    using Data = std::tuple<std::vector<typename Loggers::Type>...>;
    template <class Logger>
    static void write(const typename Logger::Type& t, Data& data)
    {
        //std::get<mio::details::IndexPosition<Logger, Types>::value>(data).push_back(t);
        std::cout << t << std::endl;
    }
};
template <class... Loggers>
struct DataWriter {
    using Data  = std::tuple<std::vector<typename Loggers::Type>...>;
    using Types = mio::Index<Loggers...>;
    template <class Logger>
    static void write(const typename Logger::Type& t, Data& data)
    {
        std::get<mio::details::IndexPosition<Logger, Types>::value>(data).push_back(t);
    }
};

template <class... Loggers>
struct ToFileWriter {
    using Data  = std::tuple<std::vector<typename Loggers::Type>...>;
    using Types = mio::Index<Loggers...>;
    template <class Logger>
    static void write(const typename Logger::Type& t, Data& data)
    {
        std::get<mio::details::IndexPosition<Logger, Types>::value>(data).push_back(t);
    }
};

void write_to_file_julia(const std::vector<int>& data)
{
    // write to file
}

template <template <class...> class Writer, class... Loggers>
class History
{
public:
    using Wri = Writer<Loggers...>;
    template <class T>
    void log(const T& t)
    {
        log_impl<T, Loggers...>(t);
    }

    const typename Wri::Data& get_log() const
    {
        return m_data;
    }

private:
    typename Wri::Data m_data;

    bool m_log_once_flag = true;

    template <class T, class logger, class... loggers>
    std::enable_if_t<std::is_base_of<LogOnce, logger>::value> log_impl(const T& t)
    {

        if (m_log_once_flag) {
            Wri::write<logger>(logger::log(t), m_data);
            m_log_once_flag = false;
        }
        log_impl<T, loggers...>(t);
    }

    template <class T, class logger, class... loggers>
    std::enable_if_t<std::is_base_of<LogAlways, logger>::value> log_impl(const T& t)
    {
        log_impl<T, loggers...>(t);
        Wri::write<logger>(logger::log(t), m_data);
    }

    template <class T>
    void log_impl(const T&)
    {
        // end of recursion
    }
};

int main()
{
    History<Writer, LogA, LogB> history;

    example ex;

    int step = 0;
    while (step < 10) {
        ex.a = step;
        ex.b = -step;

        history.log(ex);

        step++;
    }
    history.get_log();
}