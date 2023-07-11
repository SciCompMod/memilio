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

namespace details
{
/*
* @brief Helper function to get the index of a Type in a pack of Types at compile time.
* @tparam T The Type that is searched for.
* @tparam Types All Types  in the pack of Types.
* This function is used to get the index of a logger in a pack of loggers, e.g. index_templ_pack<Logger, Loggers...> gets the index of Logger in the pack Loggers.
* Only for use in a Data Writer, not at runtime.
*/
template <typename T, typename U = void, typename... Types>
constexpr size_t index_templ_pack()
{
    return std::is_same<T, U>::value ? 0 : 1 + index_templ_pack<T, Types...>();
}
} // namespace details

struct LogOnce {
};

struct LogAlways {
};

/*
* @brief This class writes data retrieved from loggers to memory. It can be used as the Writer template parameter for the History class.
* @tparam Loggers The loggers that are used to log data.
*/
template <class... Loggers>
struct DataWriterToMemory {
    using Data = std::tuple<std::vector<typename Loggers::Type>...>;
    template <class Logger>
    static void write_this(const typename Logger::Type& t, Data& data)
    {
        std::get<details::index_templ_pack<Logger, Loggers...>()>(data).push_back(t);
    }
};

/* 
* @brief History class that handles writers and Loggers.
* The class provides a log(T t) function to add the current record.
* It provides a get_log() function to access the record.
* It uses Loggers to retrieve data, and Writer to record it.
* A Logger has a type "Type", a function "static Type log(const T&)" and is derived from either LogOnce or LogAlways.
* LogOnce is only passed to Writer on the first call to History::log, LogAlways on all calls.
* The Writer defines the type "Data" of the record, and defines with "static void write_this(const Logger::Type&, Data&)" how log values are added to it.
* @tparam Writer The writer that is used to handle the data, e.g. store it into an array.
* @tparam Loggers The loggers that are used to log data.
*/

template <template <class...> class Writer, class... Loggers>
class History
{
public:
    using WriteWrapper = Writer<Loggers...>;

    template <class T>
    void log(const T& t)
    {
        log_impl<T, Loggers...>(t);
        logged = true;
    }

    /**
     * @brief Get the data object.
     * 
     * @return const WriteWrapper::Data& 
     */
    const typename WriteWrapper::Data& get_log() const
    {
        return m_data;
    }

    template <class Logger>
    const std::vector<typename Logger::Type>& get_log()
    {
        return std::get<details::index_templ_pack<Logger, Loggers...>()>(m_data);
    }

private:
    typename WriteWrapper::Data m_data;

    bool logged = false;

    template <class T, class logger, class... loggers>
    std::enable_if_t<std::is_base_of<LogOnce, logger>::value> log_impl(const T& t)
    {
        if (!logged) {
            WriteWrapper::template write_this<logger>(logger::log(t), m_data);
        }
        log_impl<T, loggers...>(t);
    }

    template <class T, class logger, class... loggers>
    std::enable_if_t<std::is_base_of<LogAlways, logger>::value> log_impl(const T& t)
    {
        log_impl<T, loggers...>(t);
        WriteWrapper::template write_this<logger>(logger::log(t), m_data);
    }

    template <class T>
    void log_impl(const T&)
    {
    }
};

template <class... Loggers>
using HistoryWithMemoryWriter = History<DataWriterToMemory, Loggers...>;

} // namespace mio
#endif //HISTORY_OBJ_H
