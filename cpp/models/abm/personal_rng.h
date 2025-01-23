/* 
* Copyright (C) 2020-2025 MEmilio
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

#ifndef MIO_ABM_PERSONAL_RNG_H
#define MIO_ABM_PERSONAL_RNG_H

#include "memilio/utils/random_number_generator.h"
#include "abm/person_id.h"
#include <cstdint>

namespace mio
{
namespace abm
{

class Person;

/**
 * Random number generator of individual persons.
 * Increments the random number generator counter of the person when used.
 * Does not store its own key or counter.
 * Instead the key needs to be provided from the outside, so that the RNG
 * for all persons share the same key.
 * The counter is taken from the person.
 * PersonalRandomNumberGenerator is cheap to construct and transparent
 * for the compiler to optimize, so we don't store the RNG persistently, only the
 * counter, so we don't need to store the key in each person. This increases
 * consistency (if the key is changed after the person is created) and
 * reduces the memory required per person.
 * @see mio::RandomNumberGeneratorBase
 */
class PersonalRandomNumberGenerator : public mio::RandomNumberGeneratorBase<PersonalRandomNumberGenerator>
{
public:
    /**
     * Creates a RandomNumberGenerator for a person.
     * @param key Key to be used by the generator.
     * @param index index of the Person.
     * @param counter Reference to the Person's RNG Counter.
     */
    PersonalRandomNumberGenerator(mio::Key<uint64_t> key, uint32_t index, mio::Counter<uint32_t>& counter);

    /**
     * Creates a RandomNumberGenerator for a person.
     * Uses the same key as another RandomNumberGenerator.
     * @param rng RandomNumberGenerator who's key will be used.
     * @param person Reference to the Person who's counter will be used.
     */
    PersonalRandomNumberGenerator(const mio::RandomNumberGenerator& rng, Person& person);

    /**
     * @return Get the key.
     */
    mio::Key<uint64_t> get_key() const
    {
        return m_key;
    }

    /**
     * @return Get the current counter.
     */
    mio::Counter<uint64_t> get_counter() const
    {
        return mio::rng_totalsequence_counter<uint64_t>(m_person_index, m_counter);
    }

    /**
     * Increment the counter.
     */
    void increment_counter()
    {
        ++m_counter;
    }

private:
    mio::Key<uint64_t> m_key; ///< Global RNG Key
    uint32_t m_person_index; ///< Index of the Person
    mio::Counter<uint32_t>& m_counter; ///< Reference to the Person's rng counter
};

} // namespace abm
} // namespace mio

#endif // MIO_ABM_PERSONAL_RNG_H
