
#ifndef ODESEIRMOBILITYIMPROVED_INFECTIONSTATE_H
#define ODESEIRMOBILITYIMPROVED_INFECTIONSTATE_H

namespace mio
{
namespace oseirmobilityimproved
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

} // namespace oseirmobilityimproved
} // namespace mio

#endif // ODESEIR_INFECTIONSTATE_H
