#ifndef EPI_ABM_AGE_H
#define EPI_ABM_AGE_H

namespace epi
{

/**
 * age groups like RKI.
 * EXPERIMENTAL; will be merged with new model framework soon.
 */
enum class AbmAgeGroup
{
    Age0to4 = 0,
    Age5to14,
    Age15to34,
    Age35to59,
    Age60to79,
    Age80plus,

    Count
};

} // namespace epi

#endif //EPI_ABM_AGE_H
