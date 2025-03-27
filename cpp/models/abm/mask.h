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

#ifndef MIO_ABM_MASK_H
#define MIO_ABM_MASK_H

#include "abm/mask_type.h"
#include "abm/time.h"
#include "memilio/io/default_serialize.h"

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
     * @param[in] t The TimePoint of the Mask's initial usage.
     */
    Mask(MaskType type, TimePoint t);

    /**
     * @brief Get the MaskType of this Mask.
     */
    MaskType get_type() const
    {
        return m_type;
    }

    /**
     * @brief Get the length of time this Mask has been used.
     * @param[in] curr_time The current TimePoint.
     */
    const TimeSpan get_time_used(TimePoint curr_time) const;

    /**
     * @brief Change the type of the Mask and reset the time it was used.
     * @param[in] new_mask_type The type of the new Mask.
     * @param[in] t The TimePoint of mask change.
     */
    void change_mask(MaskType new_mask_type, TimePoint t);

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("Mask").add("mask_type", m_type).add("time_first_used", m_time_first_usage);
    }

private:
    MaskType m_type; ///< Type of the Mask.
    TimePoint m_time_first_usage; ///< TimePoint of the Mask's initial usage.
};
} // namespace abm

/// @brief Creates an instance of abm::Mask for default serialization.
template <>
struct DefaultFactory<abm::Mask> {
    static abm::Mask create()
    {
        return abm::Mask(abm::MaskType::Count, abm::TimePoint());
    }
};

} // namespace mio

#endif
