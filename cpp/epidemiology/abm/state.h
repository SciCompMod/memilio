#ifndef EPI_ABM_STATE_H
#define EPI_ABM_STATE_H

#include <cstddef>

namespace epi
{

/** 
 * infection state in ABM.
 * can be used as 0-based index
 */
enum class InfectionState : std::size_t
{
    //TODO: More states
    Susceptible = 0,
    Exposed,
    Infected,
    Recovered,

    Count //last!!
};

}

#endif