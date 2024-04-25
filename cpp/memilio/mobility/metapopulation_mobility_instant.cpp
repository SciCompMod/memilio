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
#include "memilio/utils/compiler_diagnostics.h"

namespace mio
{
void MigrationEdge::condense_m_mobility(const double t)
{
    const size_t save_indices_size = this->m_save_indices.size();
    if (save_indices_size > 0) {

        const auto& last_value           = m_migrated.get_last_value();
        Eigen::VectorXd condensed_values = Eigen::VectorXd::Zero(save_indices_size + 1);

        // sum up the values of m_save_indices for each group (e.g. Age groups)
        std::transform(this->m_save_indices.begin(), this->m_save_indices.end(), condensed_values.data(),
                       [&last_value](const auto& indices) {
                           return std::accumulate(indices.begin(), indices.end(), 0.0,
                                                  [&last_value](double sum, auto i) {
                                                      return sum + last_value[i];
                                                  });
                       });

        // the last value is the sum of commuters
        condensed_values[save_indices_size] = m_migrated.get_last_value().sum();

        // Move the condensed values to the m_mobility_results time series
        m_mobility_results.add_time_point(t, std::move(condensed_values));
    }
}

} // namespace mio
