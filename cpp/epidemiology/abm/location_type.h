#ifndef EPI_ABM_LOCATION_TYPE_H
#define EPI_ABM_LOCATION_TYPE_H

#include <cstdint>
#include <limits>

namespace epi
{

/**
 * type of a location.
 */
enum class LocationType : std::uint32_t
{
    Home = 0,
    School,
    Work,
    SocialEvent, // TODO: differentiate different kinds
    BasicsShop, // groceries and other necessities
    Hospital,
    ICU,

    Count //last!
};

static constexpr uint32_t INVALID_LOCATION_INDEX = std::numeric_limits<uint32_t>::max();

}

#endif