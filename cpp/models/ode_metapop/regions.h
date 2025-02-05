
#ifndef ODESEIRMOBILITYIMPROVED_REGIONS_H
#define ODESEIRMOBILITYIMPROVED_REGIONS_H

#include "memilio/utils/index.h"

namespace mio
{
namespace oseirmobilityimproved
{

/**
 * @brief The Region struct is used as a dynamically
 * sized tag for all region dependent categories
 */
struct Region : public Index<Region> {
    Region(size_t val)
        : Index<Region>(val)
    {
    }
};

} // namespace oseirmobilityimproved
} // namespace mio

#endif
