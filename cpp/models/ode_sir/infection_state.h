
#ifndef SIR_INFECTIONSTATE_H
#define SIR_INFECTIONSTATE_H

namespace mio
{
namespace osir
{

/**
     * @brief The InfectionState enum describes the possible
     * categories for the infectious state of persons
     */
enum class InfectionState
{
    Susceptible,
    Infected,
    Removed,
    Count
};

} // namespace osir
} // namespace mio

#endif // SIR_INFECTIONSTATE_H
