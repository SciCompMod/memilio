/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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

#ifndef LCTSECIR_INFECTIONSTATE_H
#define LCTSECIR_INFECTIONSTATE_H

#include "memilio/utils/logging.h"
#include <vector>

namespace mio
{
namespace lsecir
{

/**
 * @brief The InfectionStateBase enum describes the possible
 * categories for the infectious state of persons
 */
enum class InfectionStateBase
{
    Susceptible        = 0,
    Exposed            = 1,
    InfectedNoSymptoms = 2,
    InfectedSymptoms   = 3,
    InfectedSevere     = 4,
    InfectedCritical   = 5,
    Recovered          = 6,
    Dead               = 7,
    Count              = 8
};

class InfectionState
{
public:
    InfectionState()
        : m_SubcompartmentNumbers(std::vector<int>((int)InfectionStateBase::Count, 1))
        , m_SubcompartmentNumbersindexfirst(std::vector<int>((int)InfectionStateBase::Count, 1))
    {
        set_compartment_index();
    }

    InfectionState(std::vector<int> SubcompartmentNumbers)
        : m_SubcompartmentNumbers(std::move(SubcompartmentNumbers))
        , m_SubcompartmentNumbersindexfirst(std::vector<int>((int)InfectionStateBase::Count, 1))
    {
        check_constraints();
        set_compartment_index();
    }

    void set_SubcompartmentNumbers(std::vector<int> SubcompartmentNumbers)
    {
        m_SubcompartmentNumbers = SubcompartmentNumbers;
        check_constraints();
        set_compartment_index();
    }

    int get_number(InfectionStateBase infectionstatebase) const
    {
        return m_SubcompartmentNumbers[(int)infectionstatebase];
    }

    int get_number(int infectionstatebaseindex) const
    {
        return m_SubcompartmentNumbers[infectionstatebaseindex];
    }

    int get_firstindex(InfectionStateBase infectionstatebase) const
    {
        return m_SubcompartmentNumbersindexfirst[(int)infectionstatebase];
    }

    int get_firstindex(int infectionstatebaseindex) const
    {
        return m_SubcompartmentNumbersindexfirst[infectionstatebaseindex];
    }

    int get_count() const
    {
        return m_Count;
    }

private:
    void check_constraints()
    {
        if (!(m_SubcompartmentNumbers.size() == (int)InfectionStateBase::Count)) {
            log_error("Vector for number of subcompartments has the wrong size.");
        }
        if (!(m_SubcompartmentNumbers[(int)InfectionStateBase::Susceptible] == 1)) {
            log_error("Susceptible compartment can not have Subcompartments.");
        }
        if (!(m_SubcompartmentNumbers[(int)InfectionStateBase::Recovered] == 1)) {
            log_error("Recovered compartment can not have Subcompartments.");
        }
        if (!(m_SubcompartmentNumbers[(int)InfectionStateBase::Dead] == 1)) {
            log_error("Dead compartment can not have Subcompartments.");
        }
        for (int i = 0; i < (int)InfectionStateBase::Count; ++i) {
            if (m_SubcompartmentNumbers[i] < 1) {
                log_error("All compartments should have at least one Subcompartment.");
            }
        }
    }
    void set_compartment_index()
    {
        int index = 0;
        for (int i = 0; i < (int)(InfectionStateBase::Count); i++) {
            m_SubcompartmentNumbersindexfirst[i] = index;
            index                                = index + m_SubcompartmentNumbers[i];
        }
        m_Count = index;
    }

    std::vector<int> m_SubcompartmentNumbers;
    std::vector<int> m_SubcompartmentNumbersindexfirst;
    int m_Count;
};

} // namespace lsecir
} // namespace mio

#endif