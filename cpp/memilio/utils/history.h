/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth
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

template <class... Types>
class Writer
{
    write(Types... args)
    {
    }
    using Data = std::tuple<Types...>;
};

template <class Writer, class... Loggers>
class History
{
public:
    //TODO: maybe change this to arglist
    template <class T>
    void log(const T& t)
    {
        log_impl(t);
        // write(Loggers::log(t)...);
    }

    TimeSeries<Writer::Data> get_data()
    {
        return m_data;
    }

private:
    Writer::Data m_data;

    template <class T, class Logger, class... Loggers>
    std::enable_if_t<std::is_derived<Logger, LogOnce>> log_impl(const T& t)
    {
        Writer::write(Logger::log(t));
        log_impl<T, Loggers...>(t);
    }

    template <class T, class Logger, class... Loggers>
    std::enable_if_t<std::is_derived<Logger, LogAlways>> log_impl(const T& t)
    {
        Writer::write(Logger::log(t));
        log_impl<T, Loggers...>(t);
    }

    template <class T>
    log_impl(const T&)
    {
        // end of recursion
    }
};
struct LogOnce {
};

struct LogLocId : public LogOnce {
    using Type = std::vector<location_id>;
    static Type log(const World& world)
    {
        [world.get_location().get_id()];
    }
};

} // namespace abm
} // namespace mio

#endif