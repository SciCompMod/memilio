#ifndef LOCATION_SPLIT_CONSTANTS_H
#define LOCATION_SPLIT_CONSTANTS_H

#include "memilio/epidemiology/age_group.h"
#include <cstddef>

/// The simulation resolves six age groups. Declarations only, see constants.cpp.
extern const size_t num_age_groups;
extern const mio::AgeGroup age_group_0_to_4;
extern const mio::AgeGroup age_group_5_to_14;
extern const mio::AgeGroup age_group_15_to_34;
extern const mio::AgeGroup age_group_35_to_59;
extern const mio::AgeGroup age_group_60_to_79;
extern const mio::AgeGroup age_group_80_plus;

#endif // LOCATION_SPLIT_CONSTANTS_H
