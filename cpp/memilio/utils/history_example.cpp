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

struct LogAlways {
};

struct LogA : LogOnce {
    using Type = int;
    static Type log(const example& ex)
    {
        return 90;
    }
};

struct LogB : LogAlways {
    using Type = int;
    static Type log(const example& ex)
    {
        return 85;
    }
};

struct LogC : LogAlways {
    using Type = std::pair<int, int>;
    static Type log(const example& ex)
    {
        return {ex.a, ex.b};
    }
};

struct LogD : LogOnce {
    using Type = int;
    static Type log(const example& ex)
    {
        return 99;
    }
};

template <class... Loggers>
struct WriterConsoleInt {
    using Data = std::tuple<std::vector<typename Loggers::Type>...>;
    template <class Logger>
    static void write_this(const typename Logger::Type& t, const Data& data)
    {
        std::cout << t << std::endl;
    }
};

template <class... Loggers>
struct WriterConsoleList {
    using Data = std::tuple<std::vector<typename Loggers::Type>...>;
    template <class Logger>
    static void write_this(const typename Logger::Type& t, const Data& data)
    {
        std::cout << t.first << " " << t.second << std::endl;
    }
};

template <typename T, typename U = void, typename... Types>
constexpr size_t indexte()
{
    return std::is_same<T, U>::value ? 0 : 1 + indexte<T, Types...>();
}

template <class... Loggers>
struct DataWriterToBuffer {
    using Data = std::tuple<std::vector<typename Loggers::Type>...>;
    template <class Logger>
    static void write_this(const typename Logger::Type& t, Data& data)
    {
        std::get<indexte<Logger, Loggers...>()>(data).push_back(t);
    }
};

void write_data_to_file(const std::vector<std::pair<int, int>>& data)
{
    for (const auto& d : data) {
        std::cout << d.first << " " << d.second << std::endl;
    }
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
    bool logged          = false;
    bool first_recursion = true;

    int index_counter = sizeof...(Loggers);

    typename Wri::Data m_data;

    template <class T, class logger, class... loggers>
    std::enable_if_t<std::is_base_of<LogOnce, logger>::value> log_impl(const T& t)
    {
        if (!logged) {
            Wri::template write_this<logger>(logger::log(t), m_data);
        }
        log_impl<T, loggers...>(t);
    }

    template <class T, class logger, class... loggers>
    std::enable_if_t<std::is_base_of<LogAlways, logger>::value> log_impl(const T& t)
    {
        log_impl<T, loggers...>(t);
        Wri::template write_this<logger>(logger::log(t), m_data);
    }

    template <class T>
    void log_impl(const T&)
    {

        if (first_recursion) { //need this to make sure every "once" logger is written just once
            logged          = true;
            first_recursion = false;
        }
        index_counter = sizeof...(Loggers) - 1;
    }
};

int main()
{
    History<WriterConsoleInt, LogA, LogB, LogD> history;
    History<WriterConsoleList, LogC> history2;
    History<DataWriterToBuffer, LogA, LogB, LogD> history3;

    example ex;

    int step = 0;
    while (step < 10) {
        ex.a = step;
        ex.b = -step;

        history.log(ex);
        history2.log(ex);
        history3.log(ex);

        step++;
    }

    auto logData = history3.get_log();
    std::cout << "EndFile and debug." << std::endl;
}