#ifndef EPI_ABM_STATE_H
#define EPI_ABM_STATE_H

#include <cstdint>

namespace epi
{

/** 
 * infection state in ABM.
 * can be used as 0-based index
 */
enum class InfectionState : std::uint32_t
{
    Susceptible = 0,
    Exposed,
    Carrier,
    Infected_Detected,
    Infected_Undetected,
    Recovered_Carrier,
    Recovered_Infected,
    Dead,

    Count //last!!
};

} // namespace epi

#endif