/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
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
#ifndef EPI_ABM_OUTPUT_RESULTS_H
#define EPI_ABM_OUTPUT_RESULTS_H

#include "abm/time.h"
#include "memilio/utils/time_series.h"
#include "abm/state.h"
#include "abm/world.h"

namespace mio
{

namespace abm
{

class OutputResults
{
public:
    OutputResults()
        : m_result(Eigen::Index(InfectionState::Count))
    {
        for (int location = (int)mio::abm::LocationType::Home; location < (int)mio::abm::LocationType::Count;
             ++location)
            m_location_result.insert(std::map<std::uint32_t, TimeSeries<double>>::value_type(
                location, TimeSeries<double>(Eigen::Index(InfectionState::Count))));
    };

    //void store_result_at(TimePoint t);
    void store_result_at(const TimePoint t, const World& world);

    void set_print_results(const bool print_results, const bool print_location_results);

    void print_result_to_file(const std::string& name_of_file);

    const char* location_type_to_string(const mio::abm::LocationType l);

private:
    TimeSeries<double> m_result;
    std::map<std::uint32_t, TimeSeries<double>> m_location_result;
    bool m_print_results          = true;
    bool m_print_location_results = true;
};
} // namespace abm
} // namespace mio

#endif // EPI_ABM_OUTPUT_RESULTS_H