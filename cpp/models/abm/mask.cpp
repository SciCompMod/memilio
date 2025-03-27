/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Carlotta Gerstein, Martin J. Kuehn
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

#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "abm/mask.h"
#include "abm/time.h"

namespace mio
{
namespace abm
{
Mask::Mask(MaskType type, TimePoint t)
    : m_type(type)
    , m_time_first_usage(t)
{
}

void Mask::change_mask(MaskType new_mask_type, TimePoint t)
{
    m_type             = new_mask_type;
    m_time_first_usage = t;
}

const TimeSpan Mask::get_time_used(TimePoint curr_time) const
{
    return curr_time - m_time_first_usage;
}

} // namespace abm
} // namespace mio
