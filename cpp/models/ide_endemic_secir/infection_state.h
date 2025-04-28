#ifndef IDE_END_SECIR_INFECTIONSTATE_H
#define IDE_END_SECIR_INFECTIONSTATE_H

namespace mio
{

namespace endisecir
{

/**
 * @brief The #InfectionState enum describes the possible
 * categories for the infectious state of persons
 */
enum class InfectionState
{
    Susceptible        = 0,
    Exposed            = 1,
    InfectedNoSymptoms = 2,
    InfectedSymptoms   = 3,
    InfectedSevere     = 4,
    InfectedCritical   = 5,
    Recovered          = 6,
    Dead               = 7,
    Count              = 8
};

/**
 * @brief The #InfectionTransition enum describes the possible
 * transitions of the infectious state of persons.
 */
enum class InfectionTransition
{
    SusceptibleToExposed                 = 0,
    ExposedToInfectedNoSymptoms          = 1,
    InfectedNoSymptomsToInfectedSymptoms = 2,
    InfectedNoSymptomsToRecovered        = 3,
    InfectedSymptomsToInfectedSevere     = 4,
    InfectedSymptomsToRecovered          = 5,
    InfectedSevereToInfectedCritical     = 6,
    InfectedSevereToRecovered            = 7,
    InfectedCriticalToDead               = 8,
    InfectedCriticalToRecovered          = 9,
    Count                                = 10
};

} // namespace endisecir
} // namespace mio

#endif //IDE_END_SECIR_INFECTIONSTATE_H
