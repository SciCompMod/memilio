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
/**
 * @brief Reduces the probability that a Person becomes infected.
 * Every Person has a Mask that reduces the probability of becoming infected when wearing this Mask.
 */
class Mask
{
public:
    /**
     * @brief Construct a new Mask of a certain type.
     * @param[in] type The type of the Mask.
     */
    Mask(MaskType type);

    /**
     * @brief Get the MaskType of this Mask.
     */
    MaskType get_type() const
    {
        return m_type;
    }

    /**
     * @brief Get the length of time this Mask has been used.
     */
    const TimeSpan& get_time_used() const
    {
        return m_time_used;
    }

    /**
     * @brief Increase the time this Mask was used by a timestep.
     * @param[in] dt The length of the timestep.
     */
    void increase_time_used(TimeSpan dt)
    {
        m_time_used += dt;
    }

    /**
     * @brief Change the type of the Mask and reset the time it was used.
     * @param[in] new_mask_type The type of the new Mask.
     */
    void change_mask(MaskType new_mask_type);

private:
    MaskType m_type; ///< Type of the Mask.
    TimeSpan m_time_used; ///< Length of time the Mask has been used.
};
} // namespace abm
} // namespace mio

#endif
