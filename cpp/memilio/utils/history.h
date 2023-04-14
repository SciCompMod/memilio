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

#ifndef MIO_ABM_HISTORY_H
#define MIO_ABM_HISTORY_H
namespace mio
{
namespace abm
{

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
        if (m_log_once_flag) { //TODO: Check if this really prints/writes all and not just the first logger
            Wri::template write_this<logger>(logger::log(t), m_data);
            m_log_once_flag = false;
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
        // end of recursion
    }
};
} // namespace abm
} // namespace mio

#endif