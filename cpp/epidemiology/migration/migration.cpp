/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "epidemiology/migration/migration.h"
#include "epidemiology/secir/seir.h"

namespace epi
{
void MigrationEdge::condense_m_migrated(double t)
{
    double infectious_over_age_groups = 0;
    double carriers_over_age_groups   = 0;
    double total_over_age_groups      = 0;
    // number of elements stored at each time point divided by secir compartments
    // equals the number of age groups
    for (auto i = 0; i < m_migrated.get_num_elements() / (int)InfectionStateV::Count; i++) {
        // sum infectious and carriers over age groups
        // m_migrated is a vector containing first all compartements of age groups 1, then age group 2, etc., hence the indexing
        infectious_over_age_groups +=
            m_migrated.get_last_value()[(int)InfectionStateV::Infected + (int)InfectionStateV::Count * (int)i];
        carriers_over_age_groups +=
            m_migrated.get_last_value()[(int)InfectionStateV::Carrier + (int)InfectionStateV::Count * (int)i];
        infectious_over_age_groups +=
            m_migrated.get_last_value()[(int)InfectionStateV::InfectedV1 + (int)InfectionStateV::Count * (int)i];
        carriers_over_age_groups +=
            m_migrated.get_last_value()[(int)InfectionStateV::CarrierV1 + (int)InfectionStateV::Count * (int)i];
        infectious_over_age_groups +=
            m_migrated.get_last_value()[(int)InfectionStateV::InfectedV2 + (int)InfectionStateV::Count * (int)i];
        carriers_over_age_groups +=
            m_migrated.get_last_value()[(int)InfectionStateV::CarrierV2 + (int)InfectionStateV::Count * (int)i];
    }
    // now get the sum over all who are allowed to travel
    for (auto i = 0; i < m_migrated.get_num_elements(); i++) {
        total_over_age_groups += m_migrated.get_last_value()[i];
    }
    // as time point t which contains now the carriers, infectious and total over age groups
    m_migrated.add_time_point(t, (TimeSeries<double>::Vector(3) << infectious_over_age_groups, carriers_over_age_groups,
                                  total_over_age_groups)
                                     .finished());
}
} // namespace epi
