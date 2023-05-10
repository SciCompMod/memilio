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
#ifndef HISTORY_OBJ_H
#define HISTORY_OBJ_H

#include <vector>
#include <tuple>
#include <iostream>

namespace mio
{

/* @brief Class to log objects to a buffer.
* This class is used to log data to a buffer, where the buffer is a tuple of vectors.
* The vectors are filled with the data of the loggers.
* For this The loggers should be structs that have a static function log() that returns the data to be logged.
* log takes as input the type of the class that is logged via the history class.
* The loggers can currently be of two types: LogOnce and LogAlways.
* LogOnce loggers are only logged once at the first call, LogAlways loggers are logged every call.
* For example this is an implementation of a logger:
* struct Logger : LogOnce{
*     using Type = int;
*     static Type log(const example& example) {
*        return example.t;
*     }
* };
* Here "example" is a class with a public int t that is logged via the history class, e.g. history.log(example).
* The History class has a get_log() function that returns the buffer.
*/

struct LogOnce {
};

struct LogAlways {
};

/*
* @brief Helper function to get the index of a Type in a pack of Types at compile time.
* @tparam T The Type that is searched for.
* @tparam Types All Types  in the pack of Types.
* This function is used to get the index of a logger in a pack of loggers, e.g. index_templ_pack<Logger, Loggers...> gets the index of Logger in the pack Loggers.
*/
template <typename T, typename U = void, typename... Types>
constexpr size_t index_templ_pack()
{
    return std::is_same<T, U>::value ? 0 : 1 + index_templ_pack<T, Types...>();
}

/*
* @brief This class is used to write data to the buffer. It is used for LogAlways loggers.
* @tparam Loggers The loggers that are used to log data.
*/
template <class... Loggers>
struct DataWriterToBuffer {
    using Data = std::tuple<std::vector<typename Loggers::Type>...>;
    template <class Logger>
    static void write_this(const typename Logger::Type& t, Data& data)
    {
        std::get<index_templ_pack<Logger, Loggers...>()>(data).push_back(t);
    }
};

/* 
* @brief Class to write data to the buffer.
* The History class is used to log data into a buffer.
* For Data only logged once the data is only logged once at the first call. For this use LogOnce Loggers.
* For Data that is logged every call use LogAlways Loggers.
* Use get_log() to get the buffer and handle the data.
* @tparam Writer The writer that is used to write data to the buffer.
* @tparam Loggers The loggers that are used to log data.
*/

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

    /**
     * @brief Get the data object.
     * 
     * @return const Wri::Data& 
     */
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

} // namespace mio
#endif //HISTORY_OBJ_H
