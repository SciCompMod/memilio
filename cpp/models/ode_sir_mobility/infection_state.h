
#ifndef ODESIRMOBILITY_INFECTIONSTATE_H
#define ODESIRMOBILITY_INFECTIONSTATE_H

namespace mio
{
namespace osirmobility
{

/**
     * @brief The InfectionState enum describes the possible
     * categories for the infectious state of persons
     */
enum class InfectionState
{
    Susceptible,
    Infected,
    Recovered,
    Count
};

} // namespace osir
} // namespace mio

#endif // ODESIR_INFECTIONSTATE_H
