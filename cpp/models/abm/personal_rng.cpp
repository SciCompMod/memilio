/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann, Khoa Nguyen, Rene Schmieding
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

#include "abm/personal_rng.h"
#include "abm/person.h"

namespace mio
{
namespace abm
{
PersonalRandomNumberGenerator::PersonalRandomNumberGenerator(mio::Key<uint64_t> key, PersonId id,
                                                             mio::Counter<uint32_t>& counter)
    : m_key(key)
    , m_person_id(id)
    , m_counter(counter)
{
}

PersonalRandomNumberGenerator::PersonalRandomNumberGenerator(const mio::RandomNumberGenerator& rng, Person& person)
    : PersonalRandomNumberGenerator(rng.get_key(), person.get_id(), person.get_rng_counter())
{
}

} // namespace abm
} // namespace mio
