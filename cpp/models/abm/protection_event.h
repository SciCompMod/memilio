/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: David Kerkmann, Khoa Nguyen
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
#ifndef MIO_ABM_VACCINE_H
#define MIO_ABM_VACCINE_H

#include "memilio/io/default_serialize.h"
#include "abm/time.h"

#include <cstdint>

namespace mio
{
namespace abm
{

/** 
 * @brief #ProtectionType in ABM.
 * can be used as 0-based index
 */
enum class ProtectionType : std::uint32_t
{
    NoProtection,
    NaturalInfection,
    GenericVaccine,
    Count //last!!
};

/**
 * @brief A tuple of #ProtectionType and #TimePoint.
 * The #TimePoint describes the time of exposure and, in case of a vaccine, the time of administration of the vaccine.
 */
struct ProtectionEvent {
    ProtectionEvent(ProtectionType exposure, TimePoint t)
        : type(exposure)
        , time(t)
    {
    }

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("ProtectionEvent").add("type", type).add("time", time);
    }

    ProtectionType type;
    TimePoint time;
};

} // namespace abm

/// @brief Creates an instance of abm::ProtectionEvent for default serialization.
template <>
struct DefaultFactory<abm::ProtectionEvent> {
    static abm::ProtectionEvent create()
    {
        return abm::ProtectionEvent(abm::ProtectionType::Count, abm::TimePoint());
    }
};

} // namespace mio

#endif
