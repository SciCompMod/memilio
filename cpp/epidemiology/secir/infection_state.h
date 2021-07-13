#ifndef INFECTIONSTATE_H
#define INFECTIONSTATE_H

namespace epi {


/**
 * @brief The InfectionState enum describes the possible
 * categories for the infectious state of persons
 */
enum class InfectionState
{
    Susceptible  = 0,
    Exposed      = 1,
    Carrier      = 2,
    Infected     = 3,
    Hospitalized = 4,
    ICU          = 5,
    Recovered    = 6,
    Dead         = 7,
    Count = 8
};

}

#endif
