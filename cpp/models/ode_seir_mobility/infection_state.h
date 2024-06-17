
#ifndef ODESEIRMOBILITY_INFECTIONSTATE_H
#define ODESEIRMOBILITY_INFECTIONSTATE_H

namespace mio
{
namespace oseirmobility
{

/**
     * @brief The InfectionState enum describes the possible
     * categories for the infectious state of persons
     */
enum class InfectionState
{
    Susceptible,
    Exposed,
    Infected,
    Recovered,
    Count
};

} // namespace oseirmobility
} // namespace mio

#endif // ODESEIR_INFECTIONSTATE_H
