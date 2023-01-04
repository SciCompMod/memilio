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

    /**
    * get the type of the mask 
    */
    MaskType get_type() const
    {
        return m_type;
    }

    /**
    * get the time the mask was used
    */
    const TimeSpan& get_time_used() const
    {
        return m_time_used;
    }

    /**
    * @brief increase the time the mask was used by a timestep
    * @param dt the length of the timestep
    */
    void increase_time_used(TimeSpan dt)
    {
        m_time_used += dt;
    }

    /**
    * @brief change the type of the mask and reset the time it was used
    * @param new_mask_type the type of the new mask
    */
    void change_mask(MaskType new_mask_type);

private:
    MaskType m_type;
    TimeSpan m_time_used;
};
} // namespace abm
} // namespace mio

#endif
