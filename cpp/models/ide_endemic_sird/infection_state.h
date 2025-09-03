#ifndef IDE_END_SIRD_INFECTIONSTATE_H
#define IDE_END_SIRD_INFECTIONSTATE_H

namespace mio
{

namespace endisird
{

/**
 * @brief The #InfectionState enum describes the possible
 * categories for the infectious state of persons
 */
enum class InfectionState
{
    Susceptible = 0,
    Infected    = 1,
    Recovered   = 2,
    Dead        = 3,
    Count       = 4
};

/**
 * @brief The #InfectionTransition enum describes the possible
 * transitions of the infectious state of persons.
 */
enum class InfectionTransition
{
    SusceptibleToInfected = 0,
    InfectedToDead        = 1,
    InfectedToRecovered   = 2,
    Count                 = 3,
};

} // namespace endisird
} // namespace mio

#endif //IDE_END_SIRD_INFECTIONSTATE_H
