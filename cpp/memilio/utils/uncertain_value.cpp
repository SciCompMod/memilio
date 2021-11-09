/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin J. Kuehn, Martin Siggel, Daniel Abele
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
#include "memilio/utils/uncertain_value.h"

namespace mio
{

void UncertainValue::set_distribution(const ParameterDistribution& dist)
{
    m_dist.reset(dist.clone());
}

observer_ptr<ParameterDistribution> UncertainValue::get_distribution()
{
    return m_dist.get();
}

observer_ptr<const ParameterDistribution> UncertainValue::get_distribution() const
{
    return m_dist.get();
}

double UncertainValue::draw_sample()
{
    if (m_dist) {
        m_value = m_dist->get_sample();
    }

    return m_value;
}

} // namespace mio
