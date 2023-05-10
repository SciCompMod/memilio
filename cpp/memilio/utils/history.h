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

// This class is used to log data to a buffer, where the buffer is a tuple of vectors.
// The vectors are filled with the data of the loggers.
// For this The loggers should be structs that have a static function log() that returns the data to be logged.
// log takes as input the type of the class that is logged via the history class.
// The loggers can currently be of two types: LogOnce and LogAlways.
// LogOnce loggers are only logged once at the first call, LogAlways loggers are logged every call.
// For example this is an implementation of a logger:
// struct Logger : LogOnce{
//     using Type = int;
//     static Type log(const example& example) {
//         return example.t;
//     }
// };
// Here example is a class with a public int t that is logged via the history class, e.g. history.log(example).
// The History class has a get_log() function that returns the buffer.

struct LogOnce {
};

struct LogAlways {
};

// This helper function returns the index of a type in a parameter pack. It is used to get the index of a logger in the buffer.
template <typename T, typename U = void, typename... Types>
constexpr size_t indexte()
{
    return std::is_same<T, U>::value ? 0 : 1 + indexte<T, Types...>();
}

// This class is used to write data to the buffer. It is used for LogAlways loggers.
template <class... Loggers>
struct DataWriterToBuffer {
    using Data = std::tuple<std::vector<typename Loggers::Type>...>;
    template <class Logger>
    static void write_this(const typename Logger::Type& t, Data& data)
    {
        std::get<indexte<Logger, Loggers...>()>(data).push_back(t);
    }
};

// This class is used to write data to the buffer. It is used for LogOnce logger
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
