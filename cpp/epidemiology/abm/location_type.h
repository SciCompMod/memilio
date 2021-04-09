#ifndef EPI_ABM_LOCATION_TYPE_H
#define EPI_ABM_LOCATION_TYPE_H

namespace epi
{

/**
 * type of a location.
 */
enum class LocationType
{
    Home,
    School,
    Work,
    SocialEvent, // TODO: differentiate different kinds
    BasicsShop, // groceries and other necessities
    Hospital,
    ICU,
};

}

#endif