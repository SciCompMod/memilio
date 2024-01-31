#ifndef MIO_ABM_PERSON_ID_H_
#define MIO_ABM_PERSON_ID_H_

#include "memilio/utils/type_safe.h"
#include <limits>

namespace mio
{
namespace abm
{

struct PersonId : mio::TypeSafe<uint32_t, PersonId>, public OperatorComparison<PersonId> {
    PersonId(uint32_t id)
        : mio::TypeSafe<uint32_t, PersonId>(id)
    {
    }
    PersonId()
        : mio::TypeSafe<uint32_t, PersonId>(std::numeric_limits<uint32_t>::max())
    {
    }

    const static PersonId invalid_id()
    {
        return PersonId();
    }
};

} // namespace abm
} // namespace mio

#endif // MIO_ABM_PERSON_ID_H_