/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele
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
#include "memilio/mobility/metapopulation_mobility_instant.h"

namespace mio
{
void MigrationEdge::condense_m_mobility(const double t, const std::vector<size_t>& indices_non_symptomatic,
                                        const std::vector<size_t>& indices_symptomatic)
{

    const auto& last_value = m_migrated.get_last_value();

    auto num_INS =
        std::accumulate(indices_non_symptomatic.begin(), indices_non_symptomatic.end(), 0., [&](auto sum, auto i) {
            return sum + last_value[i];
        });

    auto num_ISy = std::accumulate(indices_symptomatic.begin(), indices_symptomatic.end(), 0., [&](auto sum, auto i) {
        return sum + last_value[i];
    });

    double total_commuters = m_migrated.get_last_value().sum();

    // as time point t which contains now the carriers, infectious and total over age groups
    m_mobility_results.add_time_point(
        t, (mio::TimeSeries<double>::Vector(3) << num_INS, num_ISy, total_commuters).finished());
}

} // namespace mio
