#ifndef AGEGROUP_H
#define AGEGROUP_H

#include "epidemiology/utils/index.h"

namespace epi {

/**
 * @brief The AgeGroup struct is used as a dynamically
 * sized tag for all age dependent categories
 */
struct AgeGroup : public Index<AgeGroup> {
    AgeGroup(size_t val) : Index<AgeGroup>(val){}
};

}

#endif
