/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Sascha Korf
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

#ifndef EPI_ABM_IMMUNITY_LEVEL_H
#define EPI_ABM_IMMUNITY_LEVEL_H

#include "abm/time.h"

namespace mio
{
namespace abm
{

class ImmunityLevel
{
public:
    void got_infected(TimePoint t)
    {
        number_of_infections += 1;
        time_of_last_infection = t;
    };

    double get_protection_factor(VirusVariant v, TimePoint t);

    double get_severity_factor(VirusVariant v, TimePoint t);

private:
    int number_of_vaccinations         = 0;
    TimePoint time_of_last_infection   = mio::abm::TimePoint(0);
    int number_of_infections           = 0;
    TimePoint time_of_last_vaccination = mio::abm::TimePoint(0);
};
} // namespace abm
} // namespace mio

#endif