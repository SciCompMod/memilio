/*
* Copyright (C) 2020-2026 MEmilio
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
#ifndef MIO_IO_HISTORY_H
#define MIO_IO_HISTORY_H

#include "memilio/utils/metaprogramming.h"
#include <vector>
#include <tuple>
#include <iostream>

namespace mio
{

/**
 * @brief This class writes data retrieved from loggers to memory. It can be used as the Writer template parameter for the History class.
 * @tparam Loggers The loggers that are used to log data.
 */
template <class... Loggers>
struct DataWriterToMemory {
    using Data = std::tuple<std::vector<typename Loggers::Type>...>;
    /**
     * @brief Adds a new record for a given log result t to data.
     * The parameter Logger is used to determine the type of the record t, as well as the data index at which
     * the record should be added to.
     * @param[in] t The result of Logger::log.
     * @param[in,out] data An instance of Data to add the record to.
     * @tparam Logger The type of the logger used to record t.
     */
    template <class Logger>
    static void add_record(const typename Logger::Type& t, Data& data)
    {
        std::get<mio::index_of_type_v<Logger, Loggers...>>(data).push_back(t);
    }
};

/**
 * @brief History class that handles writers and loggers.
 * History provides a function "log" to add a new record and a function "get_log" to access all records.
 *
 * The History class uses Loggers to retrieve data from a given input, and a Writer to record this data.
 * A Logger is a struct with a type `Type` and functions `Type log(const T&)` and `bool should_log(const T&)`.
 * All Loggers must be unique types and default construcible/destructible. Their member "should_log" indicates whether
 * to log, while "Type" and "log" determine what is logged. The input for "should_log" and "log" is the same input
 * of type T that is given to "History::log". (Note: T does not have to be a template for a Logger implementation.)
 * The Writer defines the type `Data` to store all records (i.e. the return values of Logger::log), and the function
 * `template <class Logger> static void add_record(const Logger::Type&, Data&)` to add a new record. "add_record" is
 * used whenever "History::log" was called and "Logger::should_log" is true.
 *
 * @tparam Writer The writer that is used to handle the data, e.g. store it into an array.
 * @tparam Loggers The loggers that are used to log data.
 */
template <template <class... /*Loggers*/> class Writer, class... Loggers>
class History
{
    static_assert(!has_duplicates_v<Loggers...>, "The Loggers used for a History must be unique.");

public:
    using WriteWrapper = Writer<Loggers...>;

    History() = default;

    History(typename WriteWrapper::Data data)
        : m_data(data)
    {
    }

    History(typename Loggers::Type... args)
        : m_data(std::tie(args...))
    {
    }

    /**
     * @brief Logs new records according to the Writer and Loggers.
     *
     * Calls the log_impl function for every Logger for Input t to record data.
     * @tparam T The type of the record.
     * @param[in] t The input to record.
     */
    template <class T>
    void log(const T& t)
    {
        (log_impl(t, std::get<index_of_type_v<Loggers, Loggers...>>(m_loggers)), ...);
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

private:
    typename WriteWrapper::Data m_data;
    std::tuple<Loggers...> m_loggers;

    /**
     * @brief Checks if the given logger should log. If so, adds a record of the log to m_data.
     * @param[in] t The argument given to History::log. Passed to Logger::should_log and Logger::log.
     * @param[in] logger A Logger instance.
     * @tparam Logger A logger from the list Loggers.
     */
    template <class T, class Logger>
    void log_impl(const T& t, Logger& logger)
    {
        if (logger.should_log(t)) {
            WriteWrapper::template add_record<Logger>(logger.log(t), m_data);
        }
    }
};

template <class... Loggers>
using HistoryWithMemoryWriter = History<DataWriterToMemory, Loggers...>;

/**
 * LogOnce and LogAlways can be used as a base class to write a logger for History. They each provide the function
 * `bool should_log(const T&)`, so that only the type `Type` and the function `Type log(const T&)` have to be
 * implemented for the derived logger, where T is some input type (the same type that is given to History::log).
 * LogOnce logs only on the first call to History::log, while LogAlways logs on every call.
 *
 * For any other logging behaviour, should_log has to be defined in the logger (no base class required).
 * @see History for a full list of requirements for a logger.
 * @{
 */
struct LogOnce {
    bool was_logged = false; ///< Remember if this Logger was logged already.

    /// @brief For any type T, returns true on the first call only, and false thereafter.
    template <class T>
    bool should_log(const T&)
    {
        return was_logged ? false : (was_logged = true);
    }
};

struct LogAlways {
    /// @brief Always returns true, for any type T.
    template <class T>
    constexpr bool should_log(const T&)
    {
        return true;
    }
};
/** @} */

} // namespace mio

#endif // MIO_IO_HISTORY_H
