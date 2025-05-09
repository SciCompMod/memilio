/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef MIO_TIMER_REGISTRATION_H
#define MIO_TIMER_REGISTRATION_H

#include "memilio/timer/basic_timer.h"

#include <list>
#include <string>

namespace mio
{
namespace timing
{

struct TimerRegistration {
    std::string name;
    std::string scope;
    BasicTimer& timer;
    int thread_id;
};

/**
 * @brief If scope is empty, returns name. Otherwise, concatenates scope, two colons "::" and name.
 * @param[in] name, scope Any strings.
 * @return The qualified name given by name and scope. 
 */
inline std::string qualified_name(const std::string& name, const std::string& scope)
{
    return scope.empty() ? name : scope + "::" + name;
}

/// @brief Struct with a virtual print method to allow exchanging how TimerRegistration%s are evaluated.
struct Printer {
    virtual void print(const std::list<TimerRegistration>&, std::ostream&) = 0;

    virtual ~Printer()
    {
    }
};

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_REGISTRATION_H
