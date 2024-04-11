
#ifndef ODESIRMOBILITY_REGIONS_H
#define ODESIRMOBILITY_REGIONS_H

#include "memilio/utils/index.h"

namespace mio
{
namespace osirmobility
{

/**
 * @brief The AgeGroup struct is used as a dynamically
 * sized tag for all age dependent categories
 */
struct Region : public Index<Region> {
    Region(size_t val)
        : Index<Region>(val)
    {
    }
};

} // namespace osirmobility
} // namespace mio

#endif
