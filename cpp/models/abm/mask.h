/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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

#ifndef EPI_ABM_MASK_H
#define EPI_ABM_MASK_H

#include "abm/mask_type.h"
#include "abm/time.h"

namespace mio
{
namespace abm
{
class Mask
{
public:
    Mask(MaskType type);

    MaskType get_type() const
    {
        return m_type;
    }

    TimeSpan get_time_used()
    {
        return m_time_used;
    }

    void increase_time_used(TimeSpan dt)
    {
        m_time_used += dt;
    }

    void change_mask(MaskType new_mask_type); // changes the type of the mask and sets time_used to 0

private:
    MaskType m_type;
    TimeSpan m_time_used;
};
} // namespace abm
} // namespace mio

#endif