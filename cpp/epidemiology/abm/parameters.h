#ifndef EPI_ABM_PARAMETERS_H
#define EPI_ABM_PARAMETERS_H

#include "epidemiology/abm/age.h"
#include "epidemiology/abm/time.h"
#include "epidemiology/utils/custom_index_array.h"
#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/parameter_set.h"
#include <limits>

namespace epi
{

struct IncubationPeriod
{
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count}, 1.);
    }
};

struct SusceptibleToExposedByCarrier
{
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count}, 1.);
    }
};

struct SusceptibleToExposedByInfected
{
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count}, 1.);
    }
};

struct CarrierToInfected
{
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count}, 1.);
    }
};

struct CarrierToRecovered
{
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count}, 1.);
    }
};

struct InfectedToRecovered
{
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count}, 1.);
    }
};

struct InfectedToDead
{
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count}, 1.);
    }
};

struct RecoveredToSusceptible
{
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count}, 1.);
    }
};

struct DetectInfection
{
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static Type get_default()
    {
        return Type({AbmAgeGroup::Count}, 0.5);
    }
};


/**
 * parameters of the infection that are the same everywhere within the world.
 */
using GlobalInfectionParameters = ParameterSet<IncubationPeriod,
                                               SusceptibleToExposedByCarrier,
                                               SusceptibleToExposedByInfected,
                                               CarrierToInfected,
                                               CarrierToRecovered,
                                               InfectedToRecovered,
                                               InfectedToDead,
                                               RecoveredToSusceptible,
                                               DetectInfection>;

struct DeathFactor
{
    using Type = double;
    static constexpr Type get_default()
    {
        return 1.;
    }
};

struct EffectiveContacts
{
    using Type = double;
    static constexpr Type get_default()
    {
        return std::numeric_limits<double>::max();
    }
};

/**
 * parameters of the infection that depend on the location.
 */
using LocalInfectionParameters = ParameterSet<DeathFactor,
                                              EffectiveContacts>;

/**
 * parameters that govern the migration between locations
 */
struct LockdownDate {
    using Type = TimePoint;
    static auto get_default()
    {
        return TimePoint(std::numeric_limits<int>::max());
    }
};
struct HospitalizationRate {
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static auto get_default()
    {
        return CustomIndexArray<double, AbmAgeGroup>(AbmAgeGroup::Count, 1.0);
    }
};
struct IcuRate {
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static auto get_default()
    {
        return CustomIndexArray<double, AbmAgeGroup>(AbmAgeGroup::Count, 1.0);
    }
};
struct SocialEventRate {
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static auto get_default()
    {
        return CustomIndexArray<double, AbmAgeGroup>(AbmAgeGroup::Count, 1.0);
    }
};
struct BasicShoppingRate {
    using Type = CustomIndexArray<double, AbmAgeGroup>;
    static auto get_default()
    {
        return CustomIndexArray<double, AbmAgeGroup>(AbmAgeGroup::Count, 1.0);
    }
};

/**
 * parameters that control the migration between locations.
 */
using AbmMigrationParameters =
    ParameterSet<LockdownDate, HospitalizationRate, IcuRate, SocialEventRate, BasicShoppingRate>;

} // namespace epi
#endif
