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
    : PersonalRandomNumberGenerator(rng.get_key(), person.get_person_id(), person.get_rng_counter())
{
}

} // namespace abm
} // namespace mio